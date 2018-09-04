/*
********************************************
**        基因功能网络数据，存储边的信息       **
**      所有实现都在头文件中没有对应的源文件    **
**           write by lzxZz 2018-09-01    **
********************************************
*/              

#ifndef __EDGE_H
#define __EDGE_H
#include <string>
using std::string;


class Edge
{
public:
    //禁止隐式转换，并且同时构造出key
    explicit Edge(string g1, string g2, double wt = 0.0)
        :gene1(g1), gene2(g2), weight(wt)
    {
        key = g1 + ":" + g2;
    }
    string get_key()
    {
        return key;
    }
    double get_weight()
    {
        return weight;
    }
    ~Edge(){}

private:
    string key;         //用于hashmap的key
    string gene1;
    string gene2;
    double weight;      //用于hashmap的value
    
};

#endif