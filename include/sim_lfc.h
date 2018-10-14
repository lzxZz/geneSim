// 计算相似度的文件

#ifndef __SIM_LFC_H
#define __SIM_LFC_H


#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <fstream>
#include "defs.h"
#include "anno.h"
#include "term.h"




namespace Calculator
{
    

    class LFCValue
    {
    private:
        // ec号数组，去除掉了无效的路径
        static std::vector<std::string>                                     ec_numbers;
        // ec号到包含的基因集合
        static std::unordered_map<std::string, std::set<std::string>>       ec_genes;
        // 基因相似度hash表
        static std::unordered_map<std::string, double>                      gene_value;
       
       
        //初始化数据，参数为基因相似度文件
        static void init_data(std::string);
        // 判断两个ec号表示的路径是否有相同的基因
        static bool is_interact_by_keys(std::string, std::string);
        // 计算两个路径之间的diff，参数为基因，两个路径的ec号
        static double get_diff_by_keys(std::string,std::string, std::string);
        // 获取基因相似度，参数为g1:g2
        static double get_gene_value_by_key(std::string);
        
    public:
        //计算LFC得分，参数分别为基因相似度文件，输出文件，输出到控制台，默认值为true输出
        //基因相似度文件要求是计算通过LFCValue::gene_pair_generator(outfile);计算出来的基因对
        //      并且格式为一行一对基因,一个相似度值,使用tab分割
        static void calculator(std::string, std::string, bool = true);
        
        
    //下面是生成要计算的基因对的函数.

    private:
        // 初始化数据，不初始化基因相似度数据，用于生成要计算的基因对
        static void init_data_not_value();
    public:
        //生成要计算的基因对,参数为要生成的文件,即要计算的基因对的路径.
        static void gene_pair_generator(std::string);
    };

}

#endif