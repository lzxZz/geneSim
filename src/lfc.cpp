#include "../include/sim.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <sys/io.h>

using std::ifstream;
using std::istringstream;
using std::ofstream;

// ec号数组，去除嗲了无效的值
vector<string> Calculator::LFCValue::ec_numbers;
// 基因对到基因相似度的hashmap
unordered_map<string, double> Calculator::LFCValue::gene_value;
// ec号到包含的基因集合的hashmap
unordered_map<string, set<string>> Calculator::LFCValue::ec_genes;

// 计算LFC，参数分别为基因相似度数据文件，输出文件，是否输出到控制台
// 其中，为了避免数据覆盖，输出文件必须不存在
void Calculator::LFCValue::calculator(string gene_file, string out_file, bool is_print)
{
    if (access(out_file.c_str(), 0) == 0)
    {
        std::cout << "输出文件已经存在，请重新输入" << std::endl;
        // return;
    }
    else
    {
        std::cout << "输出文件不存在，正常工作" << std::endl;
    }
    init_data(gene_file);

    ofstream out;
    out.open(out_file);
    assert(out.is_open());

    for (auto ei : ec_numbers)
    {
        double lfc = 0;
        for (auto ej : ec_numbers)
        {
            //两个路径有相同的基因则跳过计算
            if (is_interact_by_keys(ei, ej))
            {
                continue;
            }

            auto iter = ec_genes.find(ei);
            assert(iter != ec_genes.end());

            double diff = 0;
            for (auto gene : iter->second)
            {
                if (gene == "")
                {
                    continue;
                }
                diff += get_diff_by_keys(gene, ei, ej);
            }

            lfc += diff / iter->second.size();
        }

        lfc /= ec_numbers.size();
        if (is_print)
        {
            std::cout << ei << "\t\t" << lfc << std::endl;
        }
        out << ei << "\t" << lfc << std::endl;
    }

    std::cout << "计算LFC" << std::endl;
}

bool Calculator::LFCValue::is_interact_by_keys(string ei, string ej)
{
    auto iter_ei = ec_genes.find(ei);
    auto iter_ej = ec_genes.find(ej);

    assert(iter_ei != ec_genes.end());
    assert(iter_ej != ec_genes.end());

    for (auto g1 : iter_ei->second)
    {
        for (auto g2 : iter_ej->second)
        {
            if (g1 == g2)
            {
                return true;
            }
        }
    }
    return false;
}

void Calculator::LFCValue::init_data(string gene_file)
{

    set<string> tmp_ec_numbers;
    unordered_map<string, set<string>> tmp_map;
    ifstream input_file;
    input_file.open("./data/ec.tab");

    assert(input_file.is_open());

    string line;
    while (getline(input_file, line))
    {
        vector<string> infos;

        boost::split(infos, line, boost::is_any_of("\t"));

        assert(infos.size() == 5);

        tmp_ec_numbers.emplace(infos[2]);

        auto iter = tmp_map.find(infos[2]);

        

            if (iter == tmp_map.end())
            {
                set<string> tmp_set;
                tmp_set.emplace(infos[3]);
                tmp_map.emplace(std::make_pair(infos[2], tmp_set));
            }
            else
            {
                iter->second.emplace(infos[3]);
            }
        
    }

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
        ec_numbers.push_back(ec);
        ec_genes.emplace(std::make_pair(iter->first, iter->second));
    }

    //
    //
    //读取基因数据
    //
    //
    input_file.close();

    input_file.open(gene_file);
    assert(input_file.is_open());

    while (getline(input_file, line))
    {
        istringstream is(line);
        string g1, g2;
        double value;
        is >> g1 >> g2 >> value;
        gene_value.emplace(std::make_pair(g1 + ":" + g2, value));
    }
}

double Calculator::LFCValue::get_diff_by_keys(string gene, string ei, string ej)
{
    double c = 1E-10;

    double top_value = 0, bottom_value = 0;
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
        string key = gene + ":" + gi;
        if (gene > gi)
        {
            key = gi + ":" + gene;
        }
        bottom_value += (1 - get_gene_value_by_key(key) + c);
    }
    bottom_value *= iter_ej->second.size();

    for (string gj : iter_ej->second)
    {
        if (gj == "")
        {
            continue;
        }
        string key = gene + ":" + gj;
        if (gene > gj)
        {
            key = gj + ":" + gene;
        }
        
        top_value += (1 - get_gene_value_by_key(key) + c);
    }
    top_value *= iter_ei->second.size();

    return log(top_value / bottom_value);
}

double Calculator::LFCValue::get_gene_value_by_key(string key)
{
    auto iter = gene_value.find(key);
    if (iter == gene_value.end())
    {
        return 0;
    }

    return iter->second;
}