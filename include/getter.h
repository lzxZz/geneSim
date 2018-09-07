#ifndef __GETTER_H
#define __GETTER_H
#include <set>
#include <string>
#include <vector>
#include "defs.h"
#include "data.h"
using std::set;
using std::vector;
using std::string;

namespace Data
{
    
    //提供一个多对象共享一个数据的类，所有的返回值全都是只读的。
    class Getter
    {
    
    public:
        //制定两个基因本体术语ID，获取其公共祖先节点的集合。
        static const set<string> get_public_ancestor_by_id(string, string);
        
        //获取指定的术语注释的基因集合
        static const set<string> get_term_node_anno_gene_set_by_id(string);

        //根据指定的两个术语，获取其注释的基因集合
        static const set<string> get_path_anno_gene_set_by_id(string, string);

        //获取术语注释的基因，包括注释在子孙节点上的注释，不包括id所指向的术语
        static const set<string> get_child_anno_gene_set_by_id(string);

        //根据给定的分支，获取分支上所有的注释基因集合的总数
        static const int get_root_node_anno_gene_count(Name_Space);

        //获取功能网络数据
        static const double get_net_value_by_key(string);
        //获取指定id的术语的所有子节点(直系)
        static const set<string> get_child_by_id(string);
        //获取指定id的所有子孙节点
        static const set<string> get_descendant_by_id(string);

        //获取指定注释指定id术语的基因集合
        static const set<string> get_anno_gene_set_by_id(string);

        //获取指定基因所注释的基因。
        static const set<string> get_anno_id_set_by_gene(string);
        
        //获取所有的生物过程的ec号
        static const vector<string> get_ec_numbers();

        //获取指定ec号的生物过程注释的基因
        static const set<string> get_ec_genes_by_number(string);

        //判断两个ec是否有交叉的基因
        static bool is_inter_act_by_ec_number(string, string);
    private:
        static Data datas;
        static bool is_init;
        static void init();
    };
}

#endif