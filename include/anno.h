/*
********************************************
**            注释数据，存储注释条目          **
**      所有实现都在头文件中没有对应的源文件    **
**           write by lzxZz 2018-09-02    **
********************************************
*/              

#ifndef __ANNO_H
#define __ANNO_H

#include "defs.h"
#include <string>
#include <sstream>
#include <set>
using std::string;
using std::set;
using std::ostringstream;
using std::endl;
class Annotation
{
public:
    Annotation(string go, string gene, string code, Name_Space ns)
        :go_id(go), gene_name(gene), evidence_code(code), name_space(ns)
    {

    }
    string debug()
    {
        ostringstream os;
        os << "goid:\t\t" << go_id << endl;
        os << "gene_name:\t" << gene_name << endl;
        os << "evidence_code:\t" << evidence_code << endl;
        os << "name_space:\t" << static_cast<char>(name_space) << endl;
        os << "synonym size:\t" << synonym_gene.size() << endl;
        
        os << "synonym list:" << endl;
        os << debug_synonym() << endl;

        return os.str();
    }
    string debug_synonym()
    {
        ostringstream os;

        for (auto item : synonym_gene)
        {
            os << "\t" << item << endl;
        }


        return os.str();
    }


    string get_go_id()
    {
        return go_id;
    }
    string get_gene_name()
    {
        return gene_name;
    }
    string get_evidence_code()
    {
        return evidence_code;
    }
    Name_Space get_name_space()
    {
        return name_space;
    }
    set<string> &get_synonym_gene()
    {
        return synonym_gene;
    }



    ~Annotation(){}
private:
    string          go_id;
    string          gene_name;
    string          evidence_code;
    Name_Space      name_space;
    set<string>     synonym_gene;
};



#endif