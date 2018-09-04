/*
********************************************
**        基因本体信息，存储本体节点信息       **
**      所有实现都在头文件中没有对应的源文件    **
**           write by lzxZz 2018-09-02    **
********************************************
*/              
#ifndef __TERM_H
#define __TERM_H

#include <string>
#include <sstream>
#include <set>
#include "defs.h"   //一些定义
using std::string;
using std::set;
using std::endl;
using std::ostringstream;

class Term{
private:
    string          id;
    
    Name_Space      name_space;
    bool            obsolete;    //术语过时与否
    string          name;
    set<string>     part_ids;       //part_of关系的父节点id集合
    set<string>     isa_ids;         //is_a关系的父节点id集合
public:
    Term(string _id, Name_Space ns, bool obs, string _name = "")
        :id(_id), name_space(ns), obsolete(obs),name(_name)
    {

    }

    string debug()
    {
        ostringstream os;

        os << "id:\t\t" << id << endl;
        os << "name:\t\t" << name << endl;
        os << "namespace:\t" << static_cast<char>(name_space) << endl;
        os << "part size:\t" << part_ids.size() << endl;
        os << "isa size:\t" << isa_ids.size() << endl;
        os << "obsolete:\t" << obsolete << endl;

        return os.str();
    }


    string get_id()
    {
        return id;
    }
    string get_name()
    {
        return name;
    }
    
    Name_Space get_name_space()
    {
        return name_space;
    }
    bool is_obsolete()
    {
        return obsolete;
    }
    set<string>& get_part_ids()
    {
        return part_ids;
    }
    set<string>& get_isa_ids()
    {
        return isa_ids;
    }

    



    ~Term(){}
    

};

#endif