#include "../include/calc.h"
#include <set>
#include <cmath>
using std::set;

double gene_sim(string gene1, string gene2)
{
    //从缓存中加载数据


}

//返回基因到基因集合最大相似度
double term_and_set_max_sim(string term, set<string> term_set, string ignore_genes = "")
{
    double max_value = 0.0;

    for (auto term2 : term_set)
    {
        double tmp_value = term_sim(term,term2,ignore_genes);

        max_value = max_value > tmp_value ? max_value : tmp_value;

    }

    return max_value;
}

double term_sim(string term1, string term2, string ingnore_genes)
{
    //获取公共祖先集合
    set<string> parent_set;
    double max_value = 0.0;
    
    //获取t1注释的基因集合
    int term1_gene_count;

    //获取t2注释的基因集合
    int term2_gene_count;

    //根节点术语注释的所有基因的集合
    int root_gene_coutn;

    //计算d
    double d;


    for (auto p : parent_set)
    {
        //计算U
        double uabp = get_u_abp(term1, term2, p);
        //计算f
        double fabp = d * d * uabp + (1 - d * d) * sqrt(term1_gene_count * term2_gene_count);
        
        //计算h
        double hab = d * d * root_gene_coutn + (1 - d * d) * term1_gene_count > term2_gene_count ? term1_gene_count : term2_gene_count ;

        //计算p注释的基因集合
        int parent_gene_count;

        //计算sim
        double tmp_sim = (2 * log(root_gene_coutn) - 2 * log(fabp)) 
                        / (2 * log(root_gene_coutn) - (log(term1_gene_count) + log(term2_gene_count) ) );

        tmp_sim *= (1 - hab * parent_gene_count / (root_gene_coutn * root_gene_coutn ) );
        
        max_value = max_value > tmp_sim ? max_value : tmp_sim;
    }



    return max_value;
}

//此函数所有的数据都可以进行缓存。
int get_u_abp(string ta, string tb, string tp)
{
    
    //获取ta，tb，tp的注释信息
    //计算ta到tp的路径
    //计算tb到tp的路径
}

vector<string> get_anno_by_id(string id)
{
    vector<string> result;
    //获取对应术语所注释的基因，递归获取所有的子节点注释
    return result;
}
