#include "../include/eval.h"
#include "../include/getter.h"


//计算出所有的生物路径的LFC得分,并打印到控制台
void evaluator()
{
    //获取所有的有效EC路径
    const vector<string> &ec_numbers = Data::Getter::get_ec_numbers();
    //const unordered_map<string, set<string>> ecs_gene = Data::Getter::get_ec_genes_number();
    int ec_count = 0;
    for (auto item : ec_numbers)
    {
        double lfc;
        //获取对应的基因集合
        set<string> ec_gene = Data::Getter::get_ec_genes_by_number(item);
        //计算所有的不相交集合

        for (auto ej : ec_numbers)
        {
            //如果不相交，则计算diff
            if (! Data::Getter::is_inter_act_by_ec_number(item,ej))
            {
                const set<string>& genes_ei = Data::Getter::get_ec_genes_by_number(item); 
                const set<string>& genes_ej = Data::Getter::get_ec_genes_by_number(item);
                lfc += get_diff(item, genes_ei, genes_ej)  
                        / genes_ei.size();
            }
        }
        lfc /= ec_count;

        std::cout << "EC 号:" << item << "\t\t" << lfc << std::endl;
    }
}

//判断ec1和ec2对应的基因集合是否有相交
bool is_interact(string ec1,string ec2)
{
    for (auto item1 : Data::Getter::get_ec_genes_by_number(ec1))
    {
        for (auto item2 : Data::Getter::get_ec_genes_by_number(ec2))
        {
            if (item1 == item2)
            {
                return true;
            }
        }
    }

    return false;
}

double get_diff(string gene, set<string> gene1_set, set<string> gene2_set)
{
    double top_value, bottom_value;
    double c = 1.0E-10;     //拉普拉斯平滑参数

    for (auto gene_i : gene1_set)
    {
        bottom_value += (1 + c - gene_sim(gene,gene_i));
    }
    bottom_value *= gene1_set.size();

    for (auto gene_j : gene2_set)
    {
        top_value += (1 + c - gene_sim(gene, gene_j));
    }
    top_value *= gene2_set.size();

    

    return log(top_value/bottom_value);
}