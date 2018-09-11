#include "../include/calc.h"
#include <set>
#include <cmath>
#include <assert.h>
#include <initializer_list>
using std::initializer_list;
using std::set;

double gene_sim(string gene1, string gene2)
{
    //从缓存中加载数据

    set<string> t1,t2;
    t1 = Data::Getter::get_anno_id_set_by_gene(gene1);
    t2 = Data::Getter::get_anno_id_set_by_gene(gene2);

    if (t1.size() == 0 || t2.size()  == 0)
    {
        return 0;
    }

    double top_value = 0;

    for (auto term : t1)
    {
        top_value += term_and_set_max_sim(term, t2, {gene1, gene2});
    }

    for (auto term : t2)
    {
        top_value += term_and_set_max_sim(term, t1, {gene1, gene2});
    }

    return top_value / (t1.size() + t2.size());
}

//返回基因到基因集合最大相似度
double term_and_set_max_sim(string term, set<string> term_set, initializer_list<string> ignore_genes)
{
    

    double max_value = 0.0;

    for (auto term2 : term_set)
    {
        double tmp_value = term_sim(term,term2,ignore_genes);

        max_value = max_value > tmp_value ? max_value : tmp_value;

    }

    return max_value;
}


double term_sim(string term1, string term2, initializer_list<string> ingnore_genes)
{
    //获取公共祖先集合
    set<string> parent_set = Data::Getter::get_public_ancestor_by_id(term1, term2);
    double max_value = 0.0;
    
    //获取t1注释的基因集合
    set<string> t1_annos = Data::Getter::get_anno_gene_set_by_id(term1);
    

    //获取t2注释的基因集合
    set<string> t2_annos = Data::Getter::get_anno_gene_set_by_id(term2);

    //根节点术语注释的所有基因的集合
    Term t1 = Data::Getter::get_node_by_id(term1);
    Term t2 = Data::Getter::get_node_by_id(term2);
    
    //非同分支的术语，返回0
    if (t1.get_name_space() != t2.get_name_space() )
    {
        return 0;
    }

    int root_gene_count = Data::Getter::get_root_node_anno_gene_count(t1.get_name_space());
    
    //计算d
    double d = get_dab(term1, term2, ingnore_genes);

    int term1_gene_count = Data::Getter::get_anno_gene_set_by_id(term1).size();
    int term2_gene_count = Data::Getter::get_anno_gene_set_by_id(term2).size();

    if (term1_gene_count == 0 || term2_gene_count == 0)
    {
        return 0;
    }

    for (auto p : parent_set)
    {
        //计算U
        double uabp = get_u_abp(term1, term2, p);
        // double uabp = 1;
        //计算f
        double fabp = d * d * uabp + (1 - d * d) * sqrt(term1_gene_count * term2_gene_count);
        
        //计算h
        double hab = d * d * root_gene_count + (1 - d * d) * term1_gene_count > term2_gene_count ? term1_gene_count : term2_gene_count ;

        //计算p注释的基因集合
        int parent_gene_count = Data::Getter::get_anno_gene_set_by_id(p).size();

        //计算sim
        double tmp_sim = (2 * log(root_gene_count) - 2 * log(fabp)) 
                        / (2 * log(root_gene_count) - (log(term1_gene_count) + log(term2_gene_count) ) );

        tmp_sim *= (1 - hab * parent_gene_count / (root_gene_count * root_gene_count ) );
        
        max_value = max_value > tmp_sim ? max_value : tmp_sim;
    }



    return max_value;
}

//计算两个术语对应的基因集合之间的功能距离
double get_dab(string term1, string term2, initializer_list<string> ingore_genes)
{
    

    set<string> gene_set1 = Data::Getter::get_anno_gene_set_by_id(term1);
    set<string> gene_set2 = Data::Getter::get_anno_gene_set_by_id(term2);

    for (auto gene : ingore_genes)
    {
        gene_set1.erase(gene);
        gene_set2.erase(gene);
    }
    double l12,l21;

    l12 = 0;
    for (auto g1 : gene_set1)
    {
        double tmp_value = 1;   //累乘运算，初始值要是1
        for (auto g2 : gene_set2)
        {
            //string key = g1 + ":" + g2;
            tmp_value *= (1- Data::Getter::get_net_value_by_keys(g1, g2));
        }
        l12 += tmp_value;
    }

    l21 = 0;

    for (auto g1 : gene_set2)
    {
        double tmp_value = 1;   //累乘运算，初始值要是1
        for (auto g2 : gene_set1)
        {
            //string key = g1 + ":" + g2;
            tmp_value *= (1- Data::Getter::get_net_value_by_keys(g1, g2));
        }
        l12 += tmp_value;
    }


    double value = 0;
    value = (l12 + l21) / ( 2 * (gene_set1.size() + gene_set2.size()) - l12 - l21);
    
    return value;
    
}

//此函数所有的数据都可以进行缓存。
int get_u_abp(string ta, string tb, string tp)
{
    int result = 0;
    //获取ta，tb，tp的注释信息

    set<string> gene_set;

    set<string> tmp_set;
    tmp_set = Data::Getter::get_anno_gene_set_by_id(ta);
    gene_set.insert(tmp_set.begin(),tmp_set.end());

    tmp_set = Data::Getter::get_anno_gene_set_by_id(tb);
    gene_set.insert(tmp_set.begin(),tmp_set.end());

    tmp_set = Data::Getter::get_anno_gene_set_by_id(tp);
    gene_set.insert(tmp_set.begin(),tmp_set.end());



    //计算ta到tp的路径
    tmp_set = Data::Getter::get_path_term_set_by_id(ta,tp);
    gene_set.insert(tmp_set.begin(),tmp_set.end());
    //计算tb到tp的路径
    tmp_set = Data::Getter::get_path_term_set_by_id(tb,tp);
    gene_set.insert(tmp_set.begin(),tmp_set.end());

    result = gene_set.size();

    return result;
}
