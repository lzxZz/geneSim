#ifndef __SIM_TERM_H
#define __SIM_TERM_H
#include <initializer_list>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include "anno.h"
#include "term.h"
#include "defs.h"

namespace Calculator{


    class TermSim
    {
    private:
        // 数组索引使用的方法和数据
        static std::unordered_map<std::string, int>       gene_index;
        static std::vector<std::vector<double>>           net_array;  //二维数组,存储网络值
        static std::set<std::string>                          gene_names;
        static std::unordered_map<std::string, std::vector<int>>      id_gene_index_anno;
        static int get_index_by_name(std::string);
        static double get_net_value_by_index(int, int);
        static std::vector<int> get_anno_gene_index_set_by_id(std::string);
        static void init_array_data();
        // 通过索引来计算dab
        static double get_dab_via_index(std::string, std::string);

        


        //hash操作所使用的方法和数据
        static std::vector<Annotation> gaf_items;
        static std::vector<std::string>      term_pair;
        static std::unordered_map<std::string, double>    net_value;
        
        static std::unordered_map<std::string, Term> id_term;
        // id和对应的注释基因的集合
        static std::unordered_map<std::string, std::set<std::string>>    id_gene_annos;
        static std::unordered_map<std::string, std::set<std::string>>   id_path_nodes;
        static std::unordered_map<std::string, std::set<std::string>>   id_ancestor;
        //分支注释基因总数
        static int bp_anno_count;
        static int mf_anno_count;
        static int cc_anno_count;
        //初始化数据,参数为网络数据文件
        static void init_data(std::string);
        //计算dab
        static double get_dab(std::string, std::string, std::initializer_list<std::string>);
        // 获取id对应的注释基因集合
        static std::set<std::string> get_anno_gene_set_by_id(std::string);
        static double get_net_value_by_keys(std::string, std::string);
        static void init_gaf_list();
        static void init_obo_list();
        static int get_uabp(std::string, std::string, std::string);
        static std::set<std::string> get_path_node_by_ids(std::string,std::string);
        static int get_root_gene_count_by_id(std::string);
        static void init_ancestor();
        static std::set<std::string> get_public_ancestor_by_id(std::string,std::string);
    
    public:


        // 计算两个术语之间的相似度,最后一个参数为要忽略的基因列表
        static double get_term_sim_by_ids(std::string,std::string,std::initializer_list<std::string>);

        

        // 计算术语对的显示度,为基因相似度的计算做准备,输入的网络数据文件,输出文件,并发数
        static void calculator(std::string, std::string, int = 2);
        static void calculator_by_matrix(std::string, std::string, string);
    };
}

#endif