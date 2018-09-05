#include "../include/getter.h"
#include "../include/data.h"
#include <set>
#include <string>
#include <vector>
#include "defs.h"
using std::set;
using std::vector;
using std::string;

//获取指定id的公共祖先节点集合
const set<string> &Data::Getter::get_public_ancestor_by_id(string term1, string term2)
{
    
    set<string> s;
    return s;
}







//获取指定的术语注释的基因集合
const set<string>& Data::Getter::get_term_node_anno_gene_set_by_id(string id)
{

}


//根据指定的两个术语，获取其注释的基因集合
const set<string>& Data::Getter::get_path_anno_gene_set_by_id(string term_child, string term_parent)
{
    return datas.get_path_anno_gene_set_by_id(term_child,term_parent);
}


//获取术语注释的基因，包括注释在子孙节点上的注释，不包括id所指向的术语
const set<string>& Data::Getter::get_child_anno_gene_set_by_id(string id)
{
    return datas.get_child_anno_gene_set_by_id(id);
}


//根据给定的分支，获取分支上所有的注释基因集合的总数
const int Data::Getter::get_root_node_anno_gene_count(Name_Space name_space)
{
    return datas.get_root_node_anno_gene_count(name_space);

}

//获取功能网络数据,key = g1 : g2
const double Data::Getter::get_net_value_by_key(string key)
{
    return datas.get_net_value_by_key(key);

}


//获取指定id的术语的所有子节点(直系)
const set<string>& Data::Getter::get_child_by_id(string id)
{
    return datas.get_child_by_id(id);

}


//获取指定id的所有子孙节点
const set<string>& Data::Getter::get_descendant_by_id(string id)
{
    return datas.get_descendant_by_id(id);
}

//获取指定注释指定id术语的基因集合
const set<string>& Data::Getter::get_anno_gene_set_by_id(string id)
{
    return datas.get_anno_gene_set_by_id(id);
}

//获取指定基因所注释的基因。
const set<string>& Data::Getter::get_anno_id_set_by_gene(string gene_name)
{
    return datas.get_anno_id_set_by_gene(gene_name);

}

        
//获取所有的生物过程的ec号
const vector<string>& Data::Getter::get_ec_numbers()
{
    return datas.get_ec_numbers();

}

//获取指定ec号的生物过程注释的基因
const set<string>& Data::Getter::get_ec_genes_by_number(string ec_number)
{
    return datas.get_gene_set_by_ec_number(ec_number);

}

//判断两个生物过程是否有相同的基因，有相同返回true
bool Data::Getter::is_inter_act_by_ec_number(string ec1, string ec2)
{
    set<string> set1 = get_ec_genes_by_number(ec1);
    set<string> set2 = get_ec_genes_by_number(ec2);

    for (auto gene: set1)
    {
        if ( set2.find(gene) != set2.end())
        {
            return true;
        }
    }

    return false;
}
