#include "../include/eval.h"
#include "../include/getter.h"


//计算出所有的生物路径的LFC得分,并打印到控制台
void evaluator()
{
    //获取所有的有效EC路径
    const vector<string> ec_numbers = Data::Getter::get_ec_numbers();
    //const unordered_map<string, set<string>> ecs_gene = Data::Getter::get_ec_genes_number();
    //int ec_count = 0;
    
    int round = 0;  //ec号索引
    
    
    for (auto item : ec_numbers)
    {
        double lfc = 0;
        //获取对应的基因集合
        set<string> ec_gene = Data::Getter::get_ec_genes_by_number(item);
        

        std::cout << "第" << ++round << "个ec号" << "共" << ec_numbers.size() << "个" << endl;

        //跳过所有的只有一个基因的生物过程
        if (ec_gene.size() < 2)
        {
            continue;
        }
        //计算所有的不相交集合

        for (auto ej : ec_numbers)
        {
            int index = 0;

            //如果不相交，则计算diff
            if (! Data::Getter::is_inter_act_by_ec_number(item,ej))
            {
                
                const set<string>& genes_ei = Data::Getter::get_ec_genes_by_number(item); 
                const set<string>& genes_ej = Data::Getter::get_ec_genes_by_number(ej);
                
                for (auto gene : ec_gene)
                {
                    double tmp_lfc = get_diff(gene, genes_ei, genes_ej) / genes_ei.size();
                    std::cout << "第" << round << "轮第" << ++index << "/" << ec_gene.size() << "次循环" << endl;
                    // std::cout << tmp_lfc << endl;
                    lfc += tmp_lfc;
                }
                
                
                
            }
        }
        lfc /= ec_gene.size();

        std::cout << "EC 号:" << item << "\t\t" << lfc << std::endl;
    }
}



double get_diff(string gene, set<string> gene1_set, set<string> gene2_set)
{
    double top_value = 0, bottom_value = 0;
    double c = 1.0E-5;     //拉普拉斯平滑参数

    for (auto gene_i : gene1_set)
    {
        bottom_value += (1  - gene_sim(gene,gene_i));
    }
    bottom_value *= gene2_set.size();

    for (auto gene_j : gene2_set)
    {
        top_value += (1  - gene_sim(gene, gene_j));
    }
    top_value *= gene1_set.size();

    top_value += c;
    bottom_value +=c;

    return log(top_value/bottom_value);
}