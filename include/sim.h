// 计算相似度的文件

#ifndef __SIM_H
#define __SIM_H


#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <initializer_list>
#include <fstream>
#include "defs.h"
#include "anno.h"
#include "term.h"

using std::set;
using std::vector;
using std::unordered_map;
using std::string;
using std::initializer_list;



namespace Calculator
{
    class TermSim
    {
    private:
        // 数组索引使用的方法和数据
        static unordered_map<string, int>       gene_index;
        static vector<vector<double>>           net_array;  //二维数组,存储网络值
        static set<string>                          gene_names;
        static unordered_map<string, vector<int>>      id_gene_index_anno;
        static int get_index_by_name(string);
        static double get_net_value_by_index(int, int);
        static vector<int> get_anno_gene_index_set_by_id(string);
        static void init_array_data();
        // 通过索引来计算dab
        static double get_dab_via_index(string, string);

        static unordered_map<string, double>    keys_term_sim;
        static void init_file();


        //hash操作所使用的方法和数据
        static vector<Annotation> gaf_items;
        static vector<string>      term_pair;
        static unordered_map<string, double>    net_value;
        
        static unordered_map<string, Term> id_term;
        // id和对应的注释基因的集合
        static unordered_map<string, set<string>>    id_gene_annos;
        static unordered_map<string, set<string>>   id_path_nodes;
        static unordered_map<string, set<string>>   id_ancestor;
        //分支注释基因总数
        static int bp_anno_count;
        static int mf_anno_count;
        static int cc_anno_count;
        //初始化数据,参数为网络数据文件
        static void init_data(string);
        //计算dab
        static double get_dab(string, string, initializer_list<string>);
        // 获取id对应的注释基因集合
        static set<string> get_anno_gene_set_by_id(string);
        static double get_net_value_by_keys(string, string);
        static void init_gaf_list();
        static void init_obo_list();
        static int get_uabp(string, string, string);
        static set<string> get_path_node_by_ids(string, string);
        static int get_root_gene_count_by_id(string);
        static void init_ancestor();
        static set<string> get_public_ancestor_by_id(string,string);
    public:


        // 计算两个术语之间的相似度,最后一个参数为要忽略的基因列表
        static double get_term_sim_by_ids(string,string,initializer_list<string>);

        // 计算两个术语之间的相似度,最后一个参数为要忽略的基因列表
        static double get_term_sim_by_ids_from_file(string,string,initializer_list<string>);


        // 计算术语对的显示度,为基因相似度的计算做准备,输入的网络数据文件,输出文件,并发数
        static void calculator(string, string, int = 2);
    };



    class GeneSim
    {
    private:
        // 要计算的基因对
        static vector<string>                           gene_pairs;
        // 基因到注释的术语id集合的hashmap
        static unordered_map<string, set<string>>       gene_ids_annos;  
        // 初始化数据,参数为要计算的基因对
        static void init_data(string);
        // 计算两个基因的相似度
        static double get_gene_similarity_by_keys(string, string);
        // 计算基因和一个基因集合之间最大的相似度
        static double get_max_term_and_set_sim(string, set<string>, initializer_list<string>);

        static void get_gene_similarity_by_keys_gen_ids(string, string, string);
        static void get_ids_result(string, set<string>, string);
    public:
        // 基因相似度计算，参数分别为，要计算的基因对（\t分割），输出文件夹，并行线程数
        static void calculator(string, string, int = 2);
        static void general(string,string,int =2);
    };

    class LFCValue
    {
    private:
        // ec号数组，去除掉了无效的路径
        static vector<string>                          ec_numbers;
        // ec号到包含的基因集合
        static unordered_map<string, set<string>>      ec_genes;
        // 基因相似度hash表
        static unordered_map<string, double>           gene_value;
        //初始化数据，参数为基因相似度文件
        static void init_data(string);
        // 判断两个ec号表示的路径是否有相同的基因
        static bool is_interact_by_keys(string, string);
        // 计算两个路径之间的diff，参数为基因，两个路径的ec号
        static double get_diff_by_keys(string,string, string);
        // 获取基因相似度，参数为g1:g2
        static double get_gene_value_by_key(string);
        // 初始化数据，不初始化基因相似度数据，用于生成要计算的基因对
        static void init_data_not_value();
    public:
        //计算LFC得分，参数分别为基因相似度文件，输出文件，输出到控制台，默认值为true输出
        static void calculator(string, string, bool = true);
        
        //生成要计算的基因对
        static void gene_pair_generator(string);
    };

}

#endif