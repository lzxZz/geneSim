#ifndef __SIM_GENE_H
#define __SIM_GENE_H

#include <initializer_list>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>

namespace Calculator{
    class GeneSim
    {
    private:
        // 要计算的基因对
        static std::vector<std::string>                                     gene_pairs;
        // 基因到注释的术语id集合的hashmap
        static std::unordered_map<std::string, std::set<std::string>>       gene_ids_annos;  
        //术语对相似度的map,两个goid使用:分割.
        static std::unordered_map<std::string, double>    keys_term_sim;

        // 初始化数据,参数为要计算的基因对
        static void init_data(std::string,std::string);

        // 计算两个基因的相似度
        static double get_gene_similarity_by_keys( std::string, std::string);
        // 计算基因和一个基因集合之间最大的相似度
        static double get_max_term_and_set_sim( std::string, std::set<std::string>, std::initializer_list<std::string>);

       
        
        
        // //初始化数据,输入参数为要计算的基因对和已经求值完毕的术语相似度.
        // static void init_file(std::string,std::string);
    public:
        // 基因相似度计算，参数分别为，要计算的基因对（\t分割），输出文件夹，并行线程数
        static void calculator(std::string, std::string,std::string, int = 2);
        
        
        // 计算两个术语之间的相似度,最后一个参数为要忽略的基因列表
        //最后一个参数不使用,用于除一法.暂时保留
        static double get_term_sim_by_ids_from_file(std::string, std::string,std::initializer_list<std::string>);
        
        
        
        //下面是生成要计算的术语对的方法
    private:
        
        //获取要计算的术语对,参数分别为两个要计算的基因,要输出到的文件.
        static void get_gene_similarity_by_keys_gen_ids(std::string, std::string, std::string);
        
        //循环写入要计算的术语到到指定的文件,参数分别为术语,术语集合,要输出的文件.
        //文件写入方式为追加写入.
        static void get_ids_result(std::string, std::set<std::string>, std::string);


    public:
        //生成要计算的基因对,参数为要计算的基因对,输出的术语对文件
        static void generate(std::string,std::string);


    };
}


# endif