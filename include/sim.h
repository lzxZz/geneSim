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
        static double get_term_sim_by_ids(string,string,initializer_list<string>);
    };



    class GeneSim
    {
    private:
        // 要计算的基因对
        static vector<string>                           gene_pairs;
        static unordered_map<string, set<string>>       gene_ids_annos;  

        static void init_data(string);
        static double get_gene_similarity_by_keys(string, string);
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
        //初始化数据
        static void init_data(string);
        static bool is_interact_by_keys(string, string);
        static double get_diff_by_keys(string,string, string);
        static double get_gene_value_by_key(string);
        static void init_data_not_value();
    public:
        //计算LFC得分，参数分别为基因相似度文件，输出文件，输出到控制台，默认值为true输出
        static void calculator(string, string, bool = true);
        
        //生成要计算的基因对
        static void gene_pair_generator(string);
    };

}

#endif