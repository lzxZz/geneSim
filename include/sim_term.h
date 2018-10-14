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
        

        


        //hash操作所使用的方法和数据
        static std::vector<Annotation> gaf_items;
        //要计算的术语对,数据文件硬编码为result/ids.result
        static std::vector<std::string>      term_pair;
        //功能网络数据hashmap
        static std::unordered_map<std::string, double>    net_value;

        //id到Term数据结构的hashmap        
        static std::unordered_map<std::string, Term> id_term;
        // id和对应的注释基因的集合
        static std::unordered_map<std::string, std::set<std::string>>    id_gene_annos;
        static std::unordered_map<std::string, std::set<std::string>>   id_path_nodes;
        static std::unordered_map<std::string, std::set<std::string>>   id_ancestor;
        
        //分支注释基因总数
        static int bp_anno_count;
        static int mf_anno_count;
        static int cc_anno_count;
        
        
        //计算dab,彭佳杰博士论文中公式2-1,最后一个参数为除一法中要去掉的基因,目前暂时忽略
        static double get_dab(std::string, std::string, std::initializer_list<std::string>);

        // 获取id对应的注释基因集合
        static std::set<std::string> get_anno_gene_set_by_id(std::string);
        
        //通过基因的名称搜索hashmap获取功能网络数据,速度较慢,建议使用索引获取.
        static double get_net_value_by_keys(std::string, std::string);
        
        //计算uabp数据,彭佳杰博士论文的2.2.2小节
        static int get_uabp(std::string, std::string, std::string);

        //获取两个节点之间所有的路径的所有节点,不包括路径的两个端点.
        static std::set<std::string> get_path_node_by_ids(std::string,std::string);
        
        //获取分支上所有注释基因的个数,用于计算信息量的分子.
        static int get_root_gene_count_by_id(std::string);
        
        //初始化节点到所有祖先节点的map数据,id_ancestor
        static void init_ancestor();
        //初始化注释数据,id_gene_annos
        static void init_gaf_list();
        //初始化本体数据
        static void init_obo_list();
        //初始化数据,参数为网络数据文件,硬编码要计算的术语对文件为result/ids.result.
        static void init_data(std::string);

        // 获取两个术语节点的CA集合
        static std::set<std::string> get_public_ancestor_by_id(std::string,std::string);
    

        // 计算两个术语之间的相似度,最后一个参数为要忽略的基因列表
        static double get_term_sim_by_ids(std::string,std::string,std::initializer_list<std::string>);

    public:


        
        

        // 计算术语对的相似度,为基因相似度的计算做准备,输入的网络数据文件,输出文件,并发数
        static void calculator(std::string, std::string, int = 2);
        
        //计算术语对的相似度,网络矩阵文件,基因索引映射数据文件,输出文件
        static void calculator_by_matrix(std::string, std::string, std::string);


    //++++++++++++++++++++++++++++++++++++++++++++++//
    //                                              //
    //             数组索引使用的数据和方法             //
    //                                              //
    //++++++++++++++++++++++++++++++++++++++++++++++//

    private:
        
        //基因名称到索引的hashmap
        static std::unordered_map<std::string, int>       gene_index;

        //二维数组,存储网络值
        static std::vector<std::vector<double>>           net_array;  

        //基因名称的集合(应该是一个临时数据)
        static std::set<std::string>                          gene_names;

        //id注释的基因索引集合
        static std::unordered_map<std::string, std::vector<int>>      id_gene_index_anno;

        //通过基因名称获取索引
        static int get_index_by_name(std::string);

        //通过索引获取网络值,基因1索引,基因2索引
        static double get_net_value_by_index(int, int);

        //获取id注释的基因索引集合,即访问id_gene_index_anno这一数据.
        static std::vector<int> get_anno_gene_index_set_by_id(std::string);

        //初始化网络数据,硬编码读取`data/net.txt`中的数据.
        static void init_array_data();

        // 通过索引来计算dab
        static double get_dab_via_index(std::string, std::string);
    };
}

#endif