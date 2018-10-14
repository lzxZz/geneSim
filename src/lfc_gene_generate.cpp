/*
本文件为LFCValue类中的生成要计算的基因对的实现部分.

与正常的计算逻辑上是无关的,因此拆分为两个文件.





*/

#include "../include/sim_lfc.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

using std::set;
using std::vector;
using std::unordered_map;
using std::string;


using std::ifstream;
using std::istringstream;
using std::ofstream;

// 初始化数据，跳过基因数据,配合生成要计算的基因对
void Calculator::LFCValue::init_data_not_value()
{
    //清空数据,防止重复调用时导致的错误.
    ec_numbers.clear();
    ec_genes.clear();
    
    
    set<string> tmp_ec_numbers; //为了方便的过滤无效的ec号(路径中只有一个基因)

    unordered_map<string, set<string>> tmp_map; //临时的ec号到注释基因集合的hashmap

    ifstream input_file;
    input_file.open("./data/ec.tab");
    assert(input_file.is_open());

    //读取ec.tab文件,抽取其中的ec号和基因名称.
    string line;
    while (getline(input_file, line))
    {
        vector<string> infos;
        boost::split(infos, line, boost::is_any_of("\t"));

        //所有的条目都是5列
        assert(infos.size() == 5);

        //添加ec号到临时集合.
        tmp_ec_numbers.emplace(infos[2]);

        //添加基因到ec号指向的hashmao中
        auto iter = tmp_map.find(infos[2]);
        if (iter == tmp_map.end())
        {
            //map中没有对应的ec号,添加对应的ec号.
            set<string> tmp_set;
            tmp_set.emplace(infos[3]);
            tmp_map.emplace(std::make_pair(infos[2], tmp_set));
        }
        else
        {
            iter->second.emplace(infos[3]);
        }
    }

    //过滤不需要的ec号,及其注释的基因集合
    for (auto ec : tmp_ec_numbers)
    {
        //跳过ec号为空，以及最后一位非数字的路径
        if (ec == "" || ec.back() == '-')
        {
            continue;
        }
        //跳过基因集合数目为1的集合
        auto iter = tmp_map.find(ec);
        if (iter->second.size() < 2)
        {
            continue;
        }
        //将需要的数据添加到ec_numbers,和ec_genes
        ec_numbers.push_back(ec);
        ec_genes.emplace(std::make_pair(iter->first, iter->second));
    }
}

void Calculator::LFCValue::gene_pair_generator(string out_file)
{
    if (access(out_file.c_str(), 0) == 0)
    {
        std::cout << "输出文件已经存在，请重新输入" << std::endl;
        return;
    }
    else
    {
        std::cout << "输出文件不存在，正常工作" << std::endl;
    }
    init_data_not_value();

    ofstream out;
    out.open(out_file);
    assert(out.is_open());

    set<string> gene_pair;

    for (auto ei : ec_numbers)
    {

        for (auto ej : ec_numbers)
        {
            //两个路径有相同的基因则跳过计算
            if (is_interact_by_keys(ei, ej))
            {
                continue;
            }

            auto iter = ec_genes.find(ei);
            assert(iter != ec_genes.end());

            for (auto gene : iter->second)
            {
                if (gene == "")
                {
                    continue;
                }

                auto iter_ei = ec_genes.find(ei);
                auto iter_ej = ec_genes.find(ej);

                assert(iter_ei != ec_genes.end());
                assert(iter_ej != ec_genes.end());

                for (auto gi : iter_ei->second)
                {
                    if (gi == "")
                    {
                        continue;
                    }

                    string key = gene + "\t" + gi;
                    if (gi < key)
                    {
                        key = gi + "\t" + gene;
                    }
                    gene_pair.emplace(key);
                }

                for (auto gj : iter_ej->second)
                {
                    if (gj == "")
                    {
                        continue;
                    }

                    string key = gene + "\t" + gj;
                    if (gj < key)
                    {
                        key = gj + "\t" + gene;
                    }
                    gene_pair.emplace(key);
                }
            }
        }
    }

    for (auto key : gene_pair)
    {
        out << key << std::endl;
    }
}