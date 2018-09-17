#include "../include/sim.h"

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


vector<string>                          Calculator::TermSim::term_pair;
vector<Annotation>                      Calculator::TermSim::gaf_items;
unordered_map<string, double>           Calculator::TermSim::net_value;
unordered_map<string, set<string>>      Calculator::TermSim::id_path_nodes;
unordered_map<string, set<string>>      Calculator::TermSim::id_gene_annos;
unordered_map<string, Term>             Calculator::TermSim::id_term;
unordered_map<string, set<string>>      Calculator::TermSim::id_ancestor;


// 分支上所属的注释基因数目
int  Calculator::TermSim::bp_anno_count;
int  Calculator::TermSim::mf_anno_count;
int  Calculator::TermSim::cc_anno_count;

// 计算两个术语的相似度,  忽略指定的基因(论文中的除一法)
double Calculator::TermSim::get_term_sim_by_ids(string term1, string term2, initializer_list<string> ignore_genes)
{
    // 函数所执行的操作说明见论文

    set<string> parent_set = get_public_ancestor_by_id(term1, term2); //获取公共祖先节点
    double max_value = 0.0;

    //判断两个术语是否是属于同一分支

    if (get_root_gene_count_by_id(term1) != get_root_gene_count_by_id(term2))
    {
        return 0;
    }
    
    int root_gene_count = get_root_gene_count_by_id(term1);//获取分支所有的基因数目



    // double d = get_dab(term1,term2,{}); //    获取dab的值
    double d = get_dab_via_index(term1, term2);     //索引获取网络数据,相较于hash获取快10倍
    
    int term1_gene_count = get_anno_gene_set_by_id(term1).size();
    int term2_gene_count = get_anno_gene_set_by_id(term2).size(); //获取两个数据注释的基因总数

    if (term1_gene_count == 0 || term2_gene_count == 0)
    {
        return 0;
    }

    for (auto p : parent_set)
    {
        int uabp  = get_uabp(term1, term2, p);    //获取uabp的注释信息

        //计算f
        double fabp = d * d * uabp + (1 - d * d) * sqrt(term1_gene_count * term2_gene_count);  
        
        //计算h
        double hab = d * d * root_gene_count + (1 - d * d) * term1_gene_count > term2_gene_count ? term1_gene_count : term2_gene_count ;

        //计算p注释的基因集合
        int parent_gene_count = get_anno_gene_set_by_id(p).size();

        //计算sim
        double tmp_sim = (2 * log(root_gene_count) - 2 * log(fabp)) 
                        / (2 * log(root_gene_count) - (log(term1_gene_count) + log(term2_gene_count) ) );

        tmp_sim *= (1 - hab * parent_gene_count / (root_gene_count * root_gene_count ) );
        
        max_value = max_value > tmp_sim ? max_value : tmp_sim;
    
    }

    return max_value;
}

// 暴露给外部的接口,指定网络数据文件,计算结果,输出到相对应的输出文件中, 默认线程数为2(暂未实现)
void Calculator::TermSim::calculator(string net_file, string out_file, int thread_count)
{
    
    //数据初始化
    init_data(net_file);



    ofstream out;
    out.open(out_file);
    assert(out.is_open());

    int i =0;
    boost::timer timer;

    for (auto pair : term_pair)
    {
        if (++i > 300)
        {
            
            break;
        }
        istringstream is(pair);
        string t1,t2;
        is >> t1 >> t2;
        t1 = "GO" + t1.substr(3,7);
        t2 = "GO" + t2.substr(3,7);
        double sim = get_term_sim_by_ids(t1,t2,{});
        cout << t1 << "\t" << t2 << "\t" << sim<< endl;
        out << t1 << "\t" << t2 << "\t" << sim<< endl;
    }
    cout << timer.elapsed() <<  timer.elapsed()/300 << endl;
    out.close();
}

int Calculator::TermSim::get_root_gene_count_by_id(string id)
{
    Name_Space ns = Name_Space::UNKNOWN;

    auto iter = id_term.find(id);
    if (iter != id_term.end())
    {
        ns = iter->second.get_name_space();
    }

    switch (ns)
    {
        case Name_Space::BP:
            return bp_anno_count;
        case Name_Space::MF:
            return mf_anno_count;
        case Name_Space::CC:
            return cc_anno_count;
        case Name_Space::UNKNOWN:
            return 0;
    }
    return 0;
}

void Calculator::TermSim::init_data(string net_file)
{
    // 读取网络数据文件
    ifstream input_net;
    input_net.open(net_file);
    assert(input_net.is_open());

    string line;

    while (getline(input_net, line))
    {
        // 使用字符串流,节省了手动类型转换
        istringstream is(line);
        string g1,g2;
        double value;
        is >> g1 >> g2 >> value;
        string key = g1 + ":" + g2;
        
        net_value.emplace(std::make_pair(key,value));

    }

    input_net.close();

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
    init_array_data();

}


set<string> Calculator::TermSim::get_anno_gene_set_by_id(string id)
{
    auto iter = id_gene_annos.find(id);
    if (iter == id_gene_annos.end())
    {
        return {};
    }
    return iter->second;
}

double Calculator::TermSim::get_net_value_by_keys(string g1, string g2)
{
    
    string key = g1 + ":" + g2;
    
    auto iter = net_value.find(key);
    if (iter != net_value.end())
    {
        return iter->second;
    }
    return 0;
}

int Calculator::TermSim::get_uabp(string ta, string tb, string tp)
{
    set<string> genes;

    set<string> tmp_set;

    tmp_set = get_anno_gene_set_by_id(ta);
    genes.insert(tmp_set.begin(), tmp_set.end());

    tmp_set = get_anno_gene_set_by_id(tb);
    genes.insert(tmp_set.begin(), tmp_set.end());

    tmp_set = get_anno_gene_set_by_id(tp);
    genes.insert(tmp_set.begin(), tmp_set.end());

    set<string> tmp_nodes;
    tmp_nodes = get_path_node_by_ids(ta, tp);
    for (auto id : tmp_nodes)
    {
        tmp_set = get_anno_gene_set_by_id(id);
        genes.insert(tmp_set.begin(), tmp_set.end());
    }

    tmp_nodes = get_path_node_by_ids(tb, tp);
    for (auto id : tmp_nodes)
    {
        tmp_set = get_anno_gene_set_by_id(id);
        genes.insert(tmp_set.begin(), tmp_set.end());
    }


    return genes.size();
}

set<string> Calculator::TermSim::get_path_node_by_ids(string term_child, string term_parent)
{
    string key = term_child + ":" + term_parent;

    auto iter = id_path_nodes.find(key);

    if (iter != id_path_nodes.end())
    {
        return iter->second;
    }

    return {};
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
        l21 += tmp_value;
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

        auto iter = id_gene_annos.find(go_id);

        if (iter == id_gene_annos.end())
        {
            set<string> tmp_genes;
            tmp_genes.emplace(gene_name);
            tmp_genes.insert(tmp_synonym.begin(), tmp_synonym.end());
            
            id_gene_annos.emplace(std::make_pair(go_id, tmp_genes));
        }
        
        else
        {
            iter->second.emplace(gene_name);
            iter->second.insert(tmp_synonym.begin(), tmp_synonym.end());
        }
        

    }
    
    
}


void Calculator::TermSim::init_obo_list()
{
  

    string obo_buf = "./data/obo.buf";

    ifstream buf_reader;
    buf_reader.open(obo_buf);

    if (buf_reader.is_open())
    {
        string line;
        while (getline(buf_reader, line))
        {
            vector<string> values;
            boost::split(values, line, boost::is_any_of("|"));
            assert(values.size() == 6);
            string id = values[0];
            string name = values[1];
            Name_Space ns;
            bool obs;
            set<string> iss;
            set<string> parts;

            if (values[2] == "B")
                ns = Name_Space::BP;
            else if (values[2] == "M")
                ns = Name_Space::MF;
            else if (values[2] == "C")
                ns = Name_Space::CC;
            else
                ns = Name_Space::UNKNOWN;

            if (values[3] == "true")
            {
                obs = true;
            }
            else
            {
                obs = false;
            }

            Term term{id, ns, obs, name};
            set<string> tmp_set;
            boost::split(tmp_set, values[4], boost::is_any_of("\t"));
            term.get_isa_ids().insert(tmp_set.begin(), tmp_set.end());
            term.get_isa_ids().erase("");

            boost::split(tmp_set, values[5], boost::is_any_of("\t"), boost::token_compress_on);
            term.get_part_ids().insert(tmp_set.begin(), tmp_set.end());
            term.get_part_ids().erase("");
            id_term.emplace(std::make_pair(id, term));

            
        }
        buf_reader.close();
    }
    else
    {

        string file = "./data/onto.obo";

        deque<Term> tmp_deque;
        ifstream input_file;
        input_file.open(file);
        if (!input_file.is_open())
        {
            //std::cout << "obo file read errer!\n can't open the file" << std::endl;
        }

        //节点属性参数定义
        string id;
        string name;
        Name_Space ns;
        bool obs;
        set<string> parts;
        set<string> isas;

        regex id_reg("id: GO:\\d{7}");
        regex name_reg("name: .+");
        regex obs_reg("is_obsolete: true");
        regex ns_reg("namespace: .+");
        regex part_reg("relationship: part_of GO:\\d{7}.+");
        regex isa_reg("is_a: GO:\\d{7}.+");

        string line;
        while (getline(input_file, line))
        {
            if (line == "[Term]")
            {
                Term tmp{id, ns, obs, name};

                tmp.get_part_ids().insert(parts.begin(), parts.end());
                tmp.get_isa_ids().insert(isas.begin(), isas.end());
                tmp_deque.emplace_back(tmp);

                //std::cout <<  tmp.debug() << std::endl;

                //清空数据
                id = "";
                name = "";
                ns = Name_Space::UNKNOWN;
                obs = false;
                parts.clear();
                isas.clear();
            }
            if (regex_match(line, id_reg))
            {
                id = "GO" + line.substr(7);
            }
            if (regex_match(line, name_reg))
            {
                name = line.substr(6);
            }
            if (regex_match(line, ns_reg))
            {
                if (line == "namespace: molecular_function")
                {
                    ns = Name_Space::MF;
                }

                else if (line == "namespace: cellular_component")
                {
                    ns = Name_Space::CC;
                }
                else if (line == "namespace: biological_process")
                {
                    ns = Name_Space::BP;
                }
                else
                {
                    ns = Name_Space::UNKNOWN;
                }
            }

            if (regex_match(line, obs_reg))
            {
                obs = true;
            }

            if (regex_match(line, part_reg))
            {
                parts.emplace("GO" + line.substr(25, 7));
            }
            if (regex_match(line, isa_reg))
            {
                isas.emplace("GO" + line.substr(9, 7));
            }
        }

        //加入最后一个节点，
        Term tmp{id, ns, obs, name};

        tmp.get_part_ids().insert(parts.begin(), parts.end());
        tmp.get_isa_ids().insert(isas.begin(), isas.end());
        tmp_deque.emplace_back(tmp);

        //移除掉第一个无效节点。
        tmp_deque.pop_front();
        for (auto term : tmp_deque)
        {
             id_term.emplace(term.get_id(), term);
        }
       

        ofstream buf_writer;
        buf_writer.open(obo_buf);
        assert(buf_writer.is_open());

        for (auto term : id_term)
        {
            
            buf_writer << term.second.get_id() << "|";
            buf_writer << term.second.get_name() << "|";
            buf_writer << static_cast<char>(term.second.get_name_space()) << "|";
            if (term.second.is_obsolete())
            {
                buf_writer << "true"
                           << "|";
            }
            else
            {
                buf_writer << "false"
                           << "|";
            }

            for (auto id : term.second.get_isa_ids())
            {
                buf_writer << id << "\t";
            }
            buf_writer << "|";
            for (auto id : term.second.get_part_ids())
            {
                buf_writer << id << "\t";
            }
            buf_writer << endl;
            // Log::log({__FILE__, __func__, "写入缓存完成"});
        }
    }
    
}


void Calculator::TermSim::init_ancestor()
{
    

    

    for (auto it : id_term)
    {

        deque<string> add_list;
        add_list.push_back(it.first);

        set<string> ancestor;
        while (add_list.size() > 0)
        {
            string id = add_list[0];
            auto term_iter = id_term.find(id);
            assert(term_iter != id_term.end());
            Term term = term_iter->second;

            ancestor.insert(term.get_isa_ids().begin(), term.get_isa_ids().end());
            ancestor.insert(term.get_part_ids().begin(), term.get_part_ids().end());

            add_list.insert(add_list.end(), term.get_isa_ids().begin(), term.get_isa_ids().end());
            add_list.insert(add_list.end(), term.get_part_ids().begin(), term.get_part_ids().end());

            add_list.pop_front();
        }
        id_ancestor.emplace(make_pair(it.first, ancestor));
    }


}


set<string> Calculator::TermSim::get_public_ancestor_by_id(string term1, string term2)
{
    

    set<string> public_ancestor;

    auto it1 = id_ancestor.find(term1);

    if (it1 == id_ancestor.end())
    {
        return {};
    }

    auto it2 = id_ancestor.find(term2);

    if (it2 == id_ancestor.end())
    {
        return {};
    }

    for (auto id : it1->second)
    {
        if (it2->second.find(id) != it2->second.end())
        {
            public_ancestor.emplace(id);
        }
    }

    return public_ancestor;
}

