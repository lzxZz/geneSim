#include "../include/sim.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cmath>
#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;

vector<string> Calculator::TermSim::term_pair;
vector<Annotation> Calculator::TermSim::gaf_items;
unordered_map<string, double>    Calculator::TermSim::net_value;

int  Calculator::TermSim::bp_anno_count;
int  Calculator::TermSim::mf_anno_count;
int  Calculator::TermSim::cc_anno_count;

double Calculator::TermSim::get_term_sim_by_ids(string term1, string term2, initializer_list<string> ignore_genes)
{
    set<string> parent_set; //获取公共祖先节点
    double max_value = 0.0;

    //判断两个术语是否是属于同一分支

    int root_gene_count = 1;//获取分支所有的基因数目

    double d = 1; //    获取dab的值

    int term1_gene_count = 1;
    int term2_gene_count = 1; //获取两个数据注释的基因总数

    if (term1_gene_count == 0 || term2_gene_count == 0)
    {
        return 0;
    }

    for (auto p : parent_set)
    {
        double uabp  = 0;    //获取uabp的注释信息

        //计算f
        double fabp = d * d * uabp + (1 - d * d) * sqrt(term1_gene_count * term2_gene_count);
        
        //计算h
        double hab = d * d * root_gene_count + (1 - d * d) * term1_gene_count > term2_gene_count ? term1_gene_count : term2_gene_count ;

        //计算p注释的基因集合
        int parent_gene_count = 0;

        //计算sim
        double tmp_sim = (2 * log(root_gene_count) - 2 * log(fabp)) 
                        / (2 * log(root_gene_count) - (log(term1_gene_count) + log(term2_gene_count) ) );

        tmp_sim *= (1 - hab * parent_gene_count / (root_gene_count * root_gene_count ) );
        
        max_value = max_value > tmp_sim ? max_value : tmp_sim;
    
    }

    return max_value;
}

void Calculator::TermSim::calculator(string net_file, string out_file, int thread_count)
{
    init_data(net_file);
    for (auto pair : term_pair)
    {
        istringstream is(pair);
        string t1,t2;
        is >> t1 >> t2;
        t1 = "GO" + t1.substr(3,7);
        t2 = "GO" + t2.substr(3,7);
        cout <<get_term_sim_by_ids(t1,t2,{}) << endl;
    }
}

void Calculator::TermSim::init_data(string net_file)
{
    ifstream input_net;
    input_net.open(net_file);
    assert(input_net.is_open());

    string line;

    while (getline(input_net, line))
    {
        istringstream is(line);
        string g1,g2;
        double value;
        is >> g1 >> g2 >> value;
        string key = g1 + ":" + g2;
        
        net_value.emplace(std::make_pair(key,value));

    }

    input_net.close();

    ifstream input_pair;
    input_pair.open("./result/ids.result");
    assert(input_pair.is_open());

    while (getline(input_pair, line))
    {
        term_pair.push_back(line);
    }

    input_pair.close();

    init_gaf_list();
    

}
set<string> Calculator::TermSim::get_anno_gene_set_by_id(string)
{
    return {};
}

double Calculator::TermSim::get_net_value_by_keys(string g1, string g2)
{
    return 0;
}

//计算两个术语对应的基因集合之间的功能距离
double Calculator::TermSim::get_dab(string term1, string term2, initializer_list<string> ingore_genes)
{
    

    set<string> gene_set1 = get_anno_gene_set_by_id(term1);
    set<string> gene_set2 = get_anno_gene_set_by_id(term2);

    for (auto gene : ingore_genes)
    {
        gene_set1.erase(gene);
        gene_set2.erase(gene);
    }
    double l12,l21;

    l12 = 0;
    for (auto g1 : gene_set1)
    {
        double tmp_value = 1;   //累乘运算，初始值要是1
        for (auto g2 : gene_set2)
        {
            //string key = g1 + ":" + g2;
            tmp_value *= (1- get_net_value_by_keys(g1, g2));
        }
        l12 += tmp_value;
    }

    l21 = 0;

    for (auto g1 : gene_set2)
    {
        double tmp_value = 1;   //累乘运算，初始值要是1
        for (auto g2 : gene_set1)
        {
            //string key = g1 + ":" + g2;
            tmp_value *= (1- get_net_value_by_keys(g1, g2));
        }
        l12 += tmp_value;
    }


    double value = 0;
    value = (l12 + l21) / ( 2 * (gene_set1.size() + gene_set2.size()) - l12 - l21);
    
    return value;
    
}


void Calculator::TermSim::init_gaf_list()
{
    

    string file = "./data/gene.gaf";

    ifstream input_file;

    input_file.open(file);

    string line;

    while (getline(input_file, line))
    {

        vector<string> infos;
        boost::split(infos, line, boost::is_any_of("\t"));

        //正常的注释文件有17列
        assert(infos.size() == 17);

        string go_id;
        string gene_name;
        string evidence_code;
        Name_Space name_space = Name_Space::UNKNOWN;
        set<string> synonym;

        go_id = "GO" + infos[4].substr(3);
        gene_name = infos[2];
        evidence_code = infos[6];

        string ns = infos[8];
        if (ns == "P")
        {
            name_space = Name_Space::BP;
            bp_anno_count++;
        }
        else if (ns == "F")
        {
            name_space = Name_Space::MF;
            mf_anno_count++;
        }
        else if (ns == "C")
        {
            name_space = Name_Space::CC;
            cc_anno_count++;
        }
        vector<string> tmp_synonym;
        boost::split(tmp_synonym, infos[10], boost::is_any_of("|"));

        Annotation tmp{go_id, gene_name, evidence_code, name_space};
        tmp.get_synonym_gene().insert(tmp_synonym.begin(), tmp_synonym.end());
        gaf_items.push_back(tmp);
    }
    
    
}