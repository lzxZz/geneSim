#include "../include/getter.h"

#include <set>
#include <string>
#include <vector>
#include "defs.h"
using std::set;
using std::vector;
using std::string;

//获取指定id的
const set<string> &Data::Getter::get_public_ancestor_by_id(string, string)
{
    set<string> s;
    return s;
}


//制定两个基因本体术语ID，获取其公共祖先节点的集合。
const set<string>& get_public_ancestor(string, string)
{

}


//获取指定的术语注释的基因集合的总数
const int get_term_node_anno_gene_count(string)
{

}


//获取指定的术语注释的基因集合
const set<string>& get_term_node_anno_gene_set(string)
{

}


//根据指定的两个术语，获取其注释的基因集合
const set<string>& get_path_anno_gene_set(string, string)
{

}


//获取术语注释的基因，包括注释在子孙节点上的注释，不包括id所指向的术语
const set<string>& get_term_and_child_anno_gene_set(string)
{


}


//根据给定的分支，获取分支上所有的注释基因集合的总数
const int get_root_node_anno_gene_count(Name_Space)
{


}

//获取功能网络数据
const double get_net_value_by_key(string)
{


}


//获取指定id的术语的所有子节点(直系)
const set<string>& get_child_by_id(string)
{


}


//获取指定id的所有子孙节点
const set<string>& get_descendant_by_id(string)
{

}

//获取指定注释指定id术语的基因集合
const set<string>& get_anno_gene_set_by_id(string)
{

}

//获取指定基因所注释的基因。
const set<string>& get_term_set_by_gene(string)
{


}

        
//获取所有的生物过程的ec号
const vector<string>& get_ec_numbers()
{


}

//获取指定ec号的生物过程注释的基因
const set<string>& get_ec_genes_by_number(string)
{


}

bool is_inter_act(string, string)
{


}