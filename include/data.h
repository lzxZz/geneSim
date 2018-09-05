#ifndef __DATA_H
#define __DATA_H
/*
    数据共享的类对象
*/
#include <unordered_map>
#include <string>
#include <vector>
#include <set>
#include "term.h"
#include "anno.h"
#include "edge.h"
#include "defs.h"
using std::string;
using std::vector;
using std::unordered_map;
using std::set;  
using namespace std;
namespace Data
{
    
    
    class Data
    {
     

    private:
        vector<::std::string> ec_numbers;
        unordered_map<string, set<string>> ecs_gene_set;
        //基因功能网络 key = g1:g2, value = weight / 10  做归一化处理
        unordered_map<string, double>                net_value;

        //术语ID到术语对象 key = id value = term
        unordered_map<string, Term>                  id_term;

        unordered_map<string, set<string>>           id_ancestor;

        unordered_map<string, set<string>>           id_child;

        unordered_map<string, set<string>>           id_descendant;
        

        //三个分支的注释基因数目
        int  anno_gene_count_bp;
        int  anno_gene_count_mf;
        int  anno_gene_count_cc;

        //注释在某一个术语上的基因集合 key = go_id  value = gene_set
        unordered_map<string, set<string>>           id_anno_gene;

        //ec号对应基因集合的hash map
        unordered_map<string, set<string>>           ecs_genes;  
        //ec号列表
        vector<string>                               ec_numbers;

        //本体图路径节点 key = ta|tp value = gene_set
        unordered_map<string, set<string> >           path_genes;



    public:
        vector<std::string> get_ec_numbers();
        set<string> get_gene_set_by_ec_number(string);
        set<string> get_anno_id_set_by_gene(string);
        set<string> get_anno_gene_set_by_id(string);
        set<string> get_descendant_by_id(string);
        set<string> get_child_by_id(string);
        set<string> get_path_anno_gene_set_by_id(string,string);
        set<string> get_child_anno_gene_set_by_id(string);
        set<string> get_term_node_anno_gene_set_by_id(string);
        set<string> get_public_ancestor_by_id(string,string);
        double get_net_value_by_key(string);
        int get_root_node_anno_gene_count(Name_Space);

    };


    
}

#endif 