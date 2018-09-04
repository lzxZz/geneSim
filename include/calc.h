#ifndef __CALC_H
#define __CALC_H

#include <string>
#include <vector>
#include <set>
#include "shared.h"
using std::set;
using std::vector;  
using std::string;


//计算两个基因之间的相似度
double gene_sim(string gene1, string gene2);


//计算术语和术语集合之间的最大相似度，忽略制定的基因，逗号分割多个基因
double term_and_set_max_sim(string term, set<string> term_set, string);


//计算术语term1,term2之间的术语相似度，忽略ignore_genes指定的基因，忽略的基因使用逗号分割。
double term_sim(string term1, string term2, string ingnore_genes);

//计算给予路径约束的注释信息
int get_u_abp(string ta, string tb, string tp);

//计算两个术语的公共祖先节点
set<string> get_public_ancestor(string term1, string term2);
#endif