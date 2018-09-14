#include "../include/sim.h"

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cassert>
#include <sstream>

using std::ifstream;
using std::ofstream;
using std::istringstream;



vector<string> Calculator::GeneSim::gene_pairs;


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
    return 0;
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



}