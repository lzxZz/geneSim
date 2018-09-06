#include "../include/calc.h"
#include <set>
#include <cmath>
#include <assert.h>
using std::set;

double gene_sim(string gene1, string gene2)
{
    //从缓存中加载数据

    return 0.0;
}

// //返回基因到基因集合最大相似度
// double term_and_set_max_sim(string term, set<string> term_set, string ignore_genes = "")
// {
//     double max_value = 0.0;

//     for (auto term2 : term_set)
//     {
//         double tmp_value = term_sim(term,term2,ignore_genes);

//         max_value = max_value > tmp_value ? max_value : tmp_value;

//     }

//     return max_value;
// }

// double term_sim(string term1, string term2, string ingnore_genes)
// {
//     //获取公共祖先集合
//     set<string> parent_set = get_public_ancestor(term1, term2);
//     double max_value = 0.0;
    
//     //获取t1注释的基因集合
//     auto t1_annos = id_anno_gene.find(term1);
//     assert(t1_annos != id_anno_gene.end());
//     int term1_gene_count = t1_annos->second.size();

//     //获取t2注释的基因集合
//     auto t2_annos = id_anno_gene.find(term2);
//     assert(t2_annos != id_anno_gene.end());
//     int term2_gene_count = t2_annos->second.size();

//     //根节点术语注释的所有基因的集合

//     assert(id_term.at(term1).get_name_space() == id_term.at(term2).get_name_space() );

//     int root_gene_count ;
//     switch(id_term.at(term1).get_name_space())
//     {
//         case Name_Space::BP:
//             root_gene_count = anno_gene_count_bp;
//             break;
//         case Name_Space::CC:
//             root_gene_count = anno_gene_count_cc;
//             break;
//         case Name_Space::MF:
//             root_gene_count = anno_gene_count_mf;
//             break;
//         default:
//             root_gene_count = 0;
//     }


//     //计算d
//     double d = get_dab(term1,term2);


//     for (auto p : parent_set)
//     {
//         //计算U
//         double uabp = get_u_abp(term1, term2, p);
//         //计算f
//         double fabp = d * d * uabp + (1 - d * d) * sqrt(term1_gene_count * term2_gene_count);
        
//         //计算h
//         double hab = d * d * root_gene_count + (1 - d * d) * term1_gene_count > term2_gene_count ? term1_gene_count : term2_gene_count ;

//         //计算p注释的基因集合
//         int parent_gene_count = get_anno_by_id(p).size();

//         //计算sim
//         double tmp_sim = (2 * log(root_gene_count) - 2 * log(fabp)) 
//                         / (2 * log(root_gene_count) - (log(term1_gene_count) + log(term2_gene_count) ) );

//         tmp_sim *= (1 - hab * parent_gene_count / (root_gene_count * root_gene_count ) );
        
//         max_value = max_value > tmp_sim ? max_value : tmp_sim;
//     }



//     return max_value;
// }

// //计算两个术语对应的基因集合之间的功能距离
// double get_dab(string term1, string term2)
// {

//     return 0.0;
// }

// //此函数所有的数据都可以进行缓存。
// int get_u_abp(string ta, string tb, string tp)
// {
//     int result = 0;
//     //获取ta，tb，tp的注释信息
//     auto it = id_anno_gene.find(ta);
//     assert(it != id_anno_gene.end());
//     assert(it->second.size() != 0);
//     result += it->second.size();

//     it = id_anno_gene.find(tb);
//     assert(it != id_anno_gene.end());
//     assert(it->second.size() != 0);
//     result += it->second.size();

//     it = id_anno_gene.find(tb);
//     assert(it != id_anno_gene.end());
//     assert(it->second.size() != 0);
//     result += it->second.size();


//     //计算ta到tp的路径
//     string key_ta_tp = ta + "|" + tp;
//     it = path_genes.find(key_ta_tp);
//     assert(it != path_genes.end());
//     result += it->second.size();

//     //计算tb到tp的路径
//     string key_tb_tp = tb + "|" + tp;
//     it = path_genes.find(key_tb_tp);
//     assert(it != path_genes.end());
//     result += it->second.size();

//     return result;
// }

// set<string> get_anno_by_id(string id)
// {
//     set<string> result;
//     //获取对应术语所注释的基因，递归获取所有的子节点注释,不包括id所指向的本体
//     auto it = id_descendant.find(id);

//     for (auto item : it->second)
//     {
//         for (auto gene : id_anno_gene.at(item))
//         {
//             result.insert(result.begin(),gene);
//         }
        
//     }
    


//     return result;
// }


// set<string> get_public_ancestor(string term1, string term2)
// {
//     set<string> ancestor;
    
//     auto it1 = id_ancestor.find(term1);
//     auto it2 = id_ancestor.find(term2);

//     assert(it1 != id_ancestor.end());
//     assert(it2 != id_ancestor.end());

//     for (auto item1 : it1->second)
//     {
//         for (auto item2 : it2->second)
//         {
//             if (item1 == item2)
//             {
//                 ancestor.insert(item1);
//             }
//         }
//     }

    

//     return ancestor;
// }