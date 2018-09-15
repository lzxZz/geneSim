// 计算相似度的文件

#ifndef __SIM_H
#define __SIM_H


#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <initializer_list>

using std::set;
using std::vector;
using std::unordered_map;
using std::string;
using std::initializer_list;



namespace Calculator
{
    class TermSim
    {
    public:
        // 计算两个术语之间的相似度,最后一个参数为要忽略的基因列表
        static double get_term_sim_by_ids(string,string,initializer_list<string>);
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
    public:
        // 基因相似度计算，参数分别为，要计算的基因对（\t分割），输出文件夹，并行线程数
        static void calculator(string, string, int = 2);
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