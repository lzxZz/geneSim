#include "../include/sim_gene.h"

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cassert>
#include <sstream>
#include <boost/algorithm/string.hpp>

using std::ifstream;
using std::istringstream;
using std::ofstream;

using std::string;
using std::vector;
using std::unordered_map;
using std::set;
using std::initializer_list;

//要计算的基因对
vector<string>                              Calculator::GeneSim::gene_pairs;
//基因注释的goid集合
unordered_map<string, set<string>>          Calculator::GeneSim::gene_ids_annos;

//计算基因相似度,参数分别为:要计算的基因对的文件,结果输出,术语相似度文件,并行线程数(暂未实现多线程)
void Calculator::GeneSim::calculator(string gene_pair_file, string out_file, string term_sim_file, int thread_count)
{
    //确定输出文件不存在,以防止错误参数对已有结果的覆盖
    if (access(out_file.c_str(), 0) == 0)
    {
        std::cout << "gene相似度输出文件已存在，请重新指定文件名" << std::endl;
        return;
    }
    //打开输出文件
    ofstream out;
    out.open(out_file);
    assert(out.is_open());

    // std::cout << "计算基因相似度" << std::endl;

    //初始化数据,输入参数为要计算的基因对数据,
    init_data(gene_pair_file,term_sim_file);

    for (auto key : gene_pairs)
    {
        istringstream is(key);
        string g1, g2;
        is >> g1 >> g2;
        double sim_value = get_gene_similarity_by_keys(g1, g2);
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

double Calculator::GeneSim::get_max_term_and_set_sim(string id, set<string> terms, initializer_list<string> ignore_genes)
{
    double max_value = 0.0;

    for (auto term : terms)
    {
        //修改为从文件读取术语相似度,速度更快
        double tmp_value = 0;
        if (id < term)
            tmp_value = get_term_sim_by_ids_from_file( id, term, ignore_genes);
        else
            tmp_value = get_term_sim_by_ids_from_file(term, id, ignore_genes);
        max_value = max_value > tmp_value ? max_value : tmp_value;
    }

    return max_value;
}

void Calculator::GeneSim::init_data(string gene_pair_file, string term_sim_file)
{

    if (access(gene_pair_file.c_str(), 0) != 0)
    {
        std::cout << "基因对文件不存在，请检查路径或者使用LFCValue::gene_pair_generator()生成" << std::endl;
        return;
    }
    std::cout << "开始初始化数据" << std::endl;

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
        infos.pop_back(); //移除掉数据库名称

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
            else
            {
                iter->second.emplace(go_id);
            }
        }
    }

    ifstream term_input(term_sim_file);
    assert(term_input.is_open());
    // line;
    while (getline(term_input, line))
    {
        istringstream is(line);
        string t1, t2;
        double value;
        is >> t1 >> t2 >> value;

        keys_term_sim.emplace(std::make_pair(t1 + ":" + t1, value));
    }


}

void Calculator::GeneSim::generate(string gene_pair_file, string ids_file)
{
    if (access(ids_file.c_str(), 0) == 0)
    {
        std::cout << "待计算术语对输出文件已存在，请重新指定文件名" << std::endl;
        return;
    }

    std::cout << "计算基因相似度" << std::endl;

    init_data(gene_pair_file,"");

    for (auto key : gene_pairs)
    {
        istringstream is(key);
        string g1, g2;
        is >> g1 >> g2;
        get_gene_similarity_by_keys_gen_ids(g1, g2, ids_file);
    }

    ifstream input;
    input.open(ids_file);
    assert(input.is_open());
    set<string> pairs;
    string line;
    while (getline(input, line))
    {
        pairs.emplace(line);
    }

    input.close();
    ofstream out;
    out.open(ids_file);
    assert(out.is_open());
    for (auto line : pairs)
    {
        out << line << std::endl;
    }
    out.close();
}

void Calculator::GeneSim::get_gene_similarity_by_keys_gen_ids(string g1, string g2, string ids_file)
{
    auto iter_g1 = gene_ids_annos.find(g1);
    auto iter_g2 = gene_ids_annos.find(g2);

    if (iter_g1 == gene_ids_annos.end() || iter_g2 == gene_ids_annos.end())
    {
        return;
    }

    assert(iter_g1 != gene_ids_annos.end());
    assert(iter_g2 != gene_ids_annos.end());

    for (auto t1 : iter_g1->second)
    {
        get_ids_result(t1, iter_g2->second, ids_file);
    }

    for (auto t2 : iter_g2->second)
    {
        get_ids_result(t2, iter_g1->second, ids_file);
    }
}

void Calculator::GeneSim::get_ids_result(string id, set<string> terms, string ids_file)
{
    ofstream out;
    out.open(ids_file, std::ios_base::app);

    for (auto term : terms)
    {
        //修改为从文件读取术语相似度,速度更快
        if (id < term)
        {
            out << id << "\t" << term << std::endl;
        }
        else
        {
            out << term << "\t" << id << std::endl;
        }
    }

    out.close();
}




unordered_map<string, double> Calculator::GeneSim::keys_term_sim;
// static bool is_init;

// void Calculator::GeneSim::init_file(string date_file)
// {
//     ifstream input(date_file);
//     assert(input.is_open());
//     string line;
//     while (getline(input, line))
//     {
//         istringstream is(line);
//         string t1, t2;
//         double value;
//         is >> t1 >> t2 >> value;

//         keys_term_sim.emplace(std::make_pair(t1 + ":" + t1, value));
//     }
//     is_init = true;
// }

double Calculator::GeneSim::get_term_sim_by_ids_from_file(string term1, string term2, initializer_list<string> ignore_genes)
{
    // if (!is_init)
    // {
    //     init_file(term_sim_file);
    // }
    term1 = "GO" + term1.substr(3, 7);
    term2 = "GO" + term2.substr(3, 7);
    auto iter = keys_term_sim.find(term1 + ":" + term2);

    if (iter != keys_term_sim.end())
    {
        return iter->second;
    }

    return 0;
}