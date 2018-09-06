//#include "../include/shared.h"



//基因功能网络 key = g1:g2, value = weight / 10  做归一化处理
unordered_map<string, double>                net_value{};

//ec号对应基因集合的hash map
unordered_map<string, set<string>>           ecs_genes{};  
//ec号列表
vector<string>                               ec_numbers{};

//术语ID到术语对象 key = id value = term
unordered_map<string, Term>                  id_term{};

//初始化网络数据
void 
init_net_value()
{
    vector<Edge> edges;
    read_net_file("./data/net.txt", edges);

    for (auto edge : edges)
    {
        net_value.emplace(make_pair(edge.get_key(), edge.get_weight()));
    }
    
}

void 
init_ec_tab()
{
    vector<string> ecs;
    read_ec_file("./data/ec.tab",ecs);


    unordered_map<string, set<string>> tmp;

    for (auto item : ecs)
    {
        vector<string> infos;
        split( infos, item, is_any_of("|"));

        ec_numbers.emplace_back(infos[0]);
        
        auto it = tmp.find(infos[0]);
        if ( it != tmp.end())
        {
            it->second.emplace(infos[1]);
        }
        else
        {

        set<string> tmp_set;
        tmp_set.emplace(infos[1]);
        tmp.emplace(make_pair(infos[0], tmp_set));
        }

    }
    for (auto it = tmp.begin(); it != tmp.end();++it)
    {
        if (it->second.size() > 1 && it->first != "")
        {
         
           ecs_genes.insert(make_pair(it->first,it->second));
        }
        
       
    }
    

}

void
init_term()
{
    
    vector<Term> terms;

    read_obo_file("./data/onto.obo",terms);
    for (auto term : terms)
    {
        id_term.emplace(make_pair(term.get_id(), term));
    }

}