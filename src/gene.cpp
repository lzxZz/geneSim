#include "../include/sim.h"

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cassert>
#include <sstream>
#include <boost/algorithm/string.hpp>

using std::ifstream;
using std::ofstream;
using std::istringstream;



vector<string> Calculator::GeneSim::gene_pairs;
unordered_map<string, set<string>> Calculator::GeneSim::gene_ids_annos;

void Calculator::GeneSim::calculator(string gene_pair_file, string out_file, int thread_count)
{
    if (access(out_file.c_str(), 0) == 0)
    {
        std::cout << "输出文件已存在，请重新指定文件名" << std::endl;
        //return;
    }
    ofstream out;
    out.open(out_file);
    assert(out.is_open());

    std::cout << "计算基因相似度" << std::endl;

    init_data(gene_pair_file);

    for (auto key : gene_pairs)
    {
        istringstream is(key);
        string g1,g2;
        is >> g1 >> g2;
        double sim_value = get_gene_similarity_by_keys(g1,g2);
        out << g1 << "\t" << g2 << "\t" << sim_value << std::endl;
    }


}

double Calculator::GeneSim::get_gene_similarity_by_keys(string g1, string g2)
{
    auto iter_g1 = gene_ids_annos.find(g1);
    auto iter_g2 = gene_ids_annos.find(g2);

    if (iter_g1 == gene_ids_annos.end() || iter_g2 == gene_ids_annos.end())
    {
        return 0;
    }

    assert(iter_g1 != gene_ids_annos.end());
    assert(iter_g2 != gene_ids_annos.end());

    double top_value = 0, bottom_value = 0;
    bottom_value = iter_g1->second.size() + iter_g2->second.size();

    for (auto t1 : iter_g1->second)
    {
        top_value += get_max_term_and_set_sim(t1, iter_g2->second, {g1, g2});
    }

    for (auto t2 : iter_g2->second)
    {
        top_value += get_max_term_and_set_sim(t2, iter_g1->second, {g1, g2});
    }
    
    
    return top_value / bottom_value;
}

double  Calculator::GeneSim::get_max_term_and_set_sim(string id, set<string> terms, initializer_list<string> ignore_genes)
{
    double max_value = 0.0;

    for (auto term : terms)
    {
        double tmp_value = Calculator::TermSim::get_term_sim_by_ids(id, term, ignore_genes);
        max_value = max_value > tmp_value ? max_value : tmp_value;
    }


    return max_value;
}


void Calculator::GeneSim::init_data(string gene_pair_file)
{
    
    if (access(gene_pair_file.c_str(), 0) != 0)
    {
        std::cout << "基因对文件不存在，请检查路径或者使用LFCValue::gene_pair_generator()生成" << std::endl;
        return;
    }
    std::cout << "开始初始化数据"  << std::endl;


    


    ifstream gene_input_file;
    gene_input_file.open(gene_pair_file);

    assert(gene_input_file.is_open());

    string line;

    while (getline(gene_input_file, line))
    {
        gene_pairs.push_back(line);
    }

    gene_input_file.close();

    //初始化基因到术语集合的注释信息

    ifstream anno_file;
    anno_file.open("./data/gene.gaf");

    assert(anno_file.is_open());

    while (getline(anno_file, line))
    {
        vector<string> infos;
        boost::split(infos, line, boost::is_any_of("\t"));

        // 跳过不标准的注释信息
        if (infos.size() != 17)
        {
            continue;
        }

        //获取条目中的goid和注释的所有基因
        
        set<string> genes; //注释的所有基因

        string gene_name = infos[2];
        string go_id = infos[4];
        //基因名称的同义词
        string syn = infos[10];

        
        boost::split(infos, syn, boost::is_any_of("|"));
        infos.pop_back();  //移除掉数据库名称

        genes.emplace(gene_name);
        for (auto gene : infos)
        {
            genes.emplace(gene);
        }

        // 将gene-goid的键值对插入hashmap中
        for (auto gene : genes)
        {
            auto iter = gene_ids_annos.find(gene);
            if (iter == gene_ids_annos.end())
            {
                set<string> id_set;
                id_set.emplace(go_id);
                gene_ids_annos.emplace(std::make_pair(gene, id_set));

            }
            else{
                iter->second.emplace(go_id);
            }
        }


    }
    

    

}