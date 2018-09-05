#ifndef __SHARED_H
#define __SHARED_H

#include <string>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include "data.h"
#include "anno.h"
#include "defs.h"
#include "edge.h"
#include "term.h"

using std::string;
using std::unordered_set;
using std::set;
using std::unordered_map;
using std::pair;
using std::make_pair;



//基因功能网络 key = g1:g2, value = weight / 10  做归一化处理
extern unordered_map<string, double>                net_value;

//术语ID到术语对象 key = id value = term
extern unordered_map<string, Term>                  id_term;

extern unordered_map<string, set<string>>           id_ancestor;

extern unordered_map<string, set<string>>           id_child;

extern unordered_map<string, set<string>>           id_descendant;

//三个分支的注释基因数目
extern int  anno_gene_count_bp;
extern int  anno_gene_count_mf;
extern int  anno_gene_count_cc;

//注释在某一个术语上的基因集合 key = go_id  value = gene_set
extern unordered_map<string, set<string>>           id_anno_gene;

//ec号对应基因集合的hash map
extern unordered_map<string, set<string>>           ecs_genes;  
//ec号列表
extern vector<string>                               ec_numbers;

//本体图路径节点 key = ta|tp value = gene_set
extern unordered_map<string, set<string>>           path_genes;

void init_net_value();
void init_ec_tab();
void init_term();


class Data
{
public:
    string name;
    //需要存储的各种数据
};

class Getter
{
public:
    string get_name()
    {
        return data.name;
    }
    //获取数据的方法

private:
    static Data data;
};

#endif