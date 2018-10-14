#include "../include/sim_term.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <boost/timer.hpp>
#include <regex>
#include <deque>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::deque;
using std::regex;


using std::string;
using std::vector;
using std::unordered_map;
using std::set;
using std::initializer_list;

// 基因名称到索引的映射
unordered_map<string, int>       Calculator::TermSim::gene_index;

//二维数组,存储网络值
vector<vector<double>>           Calculator::TermSim::net_array;  
// 基因名称列表,set自动去重复
set<string>                          Calculator::TermSim::gene_names;

// 获取术语注释的基因集合的索引,由于很多的基因实际上并不在网络数据之中,因此使用vector来实现,而不是set,set会影响最终结果
//  不存在的基因索引使用-1表示
unordered_map<string, vector<int>>      Calculator::TermSim::id_gene_index_anno;


void Calculator::TermSim::calculator_by_matrix(string matrix_file, string map_file, string out_file)
{

    //确定输出文件不存在,以防止错误参数对已有结果的覆盖
    if (access(out_file.c_str(), 0) == 0)
    {
        std::cout << "术语相似度-相似度输出文件已存在，请重新指定文件名" << std::endl;
        return;
    }


    //首先完成基因名称到索引的映射.
    ifstream input_map;
    input_map.open(map_file);
    assert(input_map.is_open());

    int index = 0;

    string line;
    while (getline(input_map, line))
    {
        gene_index.emplace(std::make_pair(line, index));
        ++index;
    }

    cout << "Gene Numbers: " << gene_index.size()  << endl;

    ifstream input_matrix;
    input_matrix.open(matrix_file);
    assert(input_matrix.is_open());
    index = 0;
    while (getline(input_matrix, line))
    {
        net_array.push_back({});
        double value;
        istringstream is(line);
        while (is >> value)
        {
            net_array.at(index).push_back(value);
        }
        ++index ;
    }



    std::cout << "矩阵初始化完成" << std::endl;

    //读取要计算的术语对
    ifstream input_pair;
    input_pair.open("./result/ids.result");
    assert(input_pair.is_open());

    while (getline(input_pair, line))
    {
        term_pair.push_back(line);
    }

    input_pair.close();

    // 初始化注释文件相关数据
    init_gaf_list();

    ifstream input_path;
    input_path.open("./data/path.buf");
    assert(input_path.is_open());

    while (getline(input_path, line))
    {
        vector<string> infos;
        boost::split(infos, line, boost::is_any_of("|"));

        assert(infos.size() == 2);

        string key = infos[0];
        string values = infos[1];

        boost::split(infos, values, boost::is_any_of("\t"));

        set<string> nodes;
        nodes.insert(infos.begin(), infos.end());

        id_path_nodes.emplace(std::make_pair(key, nodes));
    }

    init_obo_list();
    init_ancestor();
    // init_array_data();

    // 构造id到基因索引的hashmap
    for (auto iter : id_gene_annos)
    {
        
        vector<int> indexs;

        for (auto gene : iter.second)
        {
            //将基因名称转换为索引
            indexs.push_back(get_index_by_name(gene));
        }

        id_gene_index_anno.emplace(std::make_pair(iter.first, indexs));
    }


    ofstream out;
    out.open(out_file);
    assert(out.is_open());

    

    int i = 0;
    boost::timer timer;

    for (auto pair : term_pair)
    {

        istringstream is(pair);
        string t1, t2;
        is >> t1 >> t2;
        t1 = "GO" + t1.substr(3, 7);
        t2 = "GO" + t2.substr(3, 7);
        double sim = get_term_sim_by_ids(t1, t2, {});
        cout << t1 << "\t" << t2 << "\t" << sim << endl;
        out << t1 << "\t" << t2 << "\t" << sim << endl;

        cout << "计算" << i++ << "对共用时" << timer.elapsed() << endl;
    }

    out.close();

}

// 初始化索引访问所需要的数据
void Calculator::TermSim::init_array_data()
{
    //开始读取网络数据
    ifstream input_net;
    input_net.open("./data/.txt");
    assert(input_net.is_open());
    
    //读取基因的名称
    string line;
    while (getline(input_net, line))
    {
        string g1,g2;
        istringstream is(line);
        is >> g1 >> g2;
        gene_names.emplace(g1);
        gene_names.emplace(g2);
    }

    //用0初始化矩阵
    for (size_t i = 0; i < gene_names.size(); i++)
    {
        vector<double> tmp_vec(gene_names.size(),0);
        
        net_array.push_back(tmp_vec);
    }

    //构造基因名称到索引的hashmap
    int index = 0;
    for (auto gene : gene_names)
    {
        gene_index.emplace(std::make_pair(gene, index));
        index++;
    }

    // 构造id到基因索引的hashmap
    for (auto iter : id_gene_annos)
    {
        
        vector<int> indexs;

        for (auto gene : iter.second)
        {
            //将基因名称转换为索引
            indexs.push_back(get_index_by_name(gene));
        }

        id_gene_index_anno.emplace(std::make_pair(iter.first, indexs));
    }

    //读取网络数据文件,覆盖掉0矩阵中相对应的值
    for (auto iter : net_value)
    {
        string key = iter.first;
        vector<string> genes;
        boost::split(genes, key, boost::is_any_of(":"));
       
        assert(genes.size() == 2);//正常的都会有两个基因名称
        int x,y;
        x = get_index_by_name(genes[0]);
        y = get_index_by_name(genes[1]);

        

        net_array.at(x).at(y) = iter.second;
    }


}

// 获取基因对应的索引,不存在的基因返回-1的索引
int Calculator::TermSim::get_index_by_name(string name)
{
    cout << gene_index.size() << endl;

    auto iter = gene_index.find(name);

    if (iter != gene_index.end())
    {
        return iter->second;
    }

    return -1;
}

// 获取术语对应的注释基因集合的索引集合
vector<int> Calculator::TermSim::get_anno_gene_index_set_by_id(string id)
{
    auto iter = id_gene_index_anno.find(id);

    if (iter != id_gene_index_anno.end())
    {
        return iter->second;
    }

    return {};
}

//  获取网络数据的值,通过索引
double Calculator::TermSim::get_net_value_by_index(int x, int y)
{
    
    if (x < 0 || y < 0 )
    {
        return 0;
    }
    if (x >= gene_index.size()|| y >= gene_index.size()) //--------------> 此处使用int和size_t比较是因为维度可能有负数,但是上面已经判断过了,是安全的
    {
        return 0;
    }
    
    
    return net_array.at(x).at(y);   
    
    
    

}


//计算两个术语对应的基因集合之间的功能距离,程序整体结构与get_adb(string, string)一致,但是在获取网络数据的时候使用索引进行访问.
double Calculator::TermSim::get_dab_via_index(string term1, string term2)
{
    //具体计算步骤参见公式201

    vector<int> gene_set1 = get_anno_gene_index_set_by_id(term1);
    vector<int> gene_set2 = get_anno_gene_index_set_by_id(term2);

   

    double l12,l21;

    l12 = 0;
    for (auto index1 : gene_set1)
    {
        double tmp_value = 1;   //累乘运算，初始值要是1
        for (auto index2 : gene_set2)
        {
            //string key = g1 + ":" + g2;
            
            tmp_value *= (1- get_net_value_by_index(index1, index2));
        }
        l12 += tmp_value;
    }

    l21 = 0;

    for (auto index1 : gene_set2)
    {
        double tmp_value = 1;   //累乘运算，初始值要是1
        for (auto index2 : gene_set1)
        {
            //string key = g1 + ":" + g2;
            tmp_value *= (1- get_net_value_by_index(index1, index2));
        }
        l21 += tmp_value;
    }


    double value = 0;
    value = (l12 + l21) / ( 2 * (gene_set1.size() + gene_set2.size()) - l12 - l21);
    
    return value;
    
}