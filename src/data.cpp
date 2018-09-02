#include "../include/data.h"
#include "../include/term.h"
#include "../include/anno.h"
#include <deque>
#include <regex>
using std::regex;
using std::deque;

void read_net_file(const string file,vector<Edge> &edges)
{
    
    ifstream input_file;
    
    input_file.open(file);
    
    string line;
    while(getline(input_file,line)){
        istringstream is(line);
        string g1,g2;
        double weight;
        is >> g1;
        is >> g2;
        is >> weight;
        edges.push_back(Edge{g1,g2,weight});
    }
}

void read_obo_file(const string file,deque<Term> &terms)
{
    ifstream input_file;
    input_file.open(file);
    if (! input_file.is_open()){
        std::cout << "obo file read errer!\n can't open the file" << std::endl;
    }

    //节点属性参数定义
    string          id;
    string          name;
    Name_Space      ns;
    bool            obs;
    set<string>     parts;
    set<string>     isas;

    regex id_reg("id: GO:\\d{7}");
    regex name_reg("name: .+");
    regex obs_reg("is_obsolete: true");
    regex ns_reg("namespace: .+");
    regex part_reg("relationship: part_of GO:\\d{7}.+");
    regex isa_reg("is_a: GO:\\d{7}.+");
    

    string line;
    while(getline(input_file,line))
    {
        if (line == "[Term]")
        {
            Term tmp{id,ns,obs,name};

            tmp.get_part_ids().insert(parts.begin(), parts.end());
            tmp.get_isa_ids().insert(isas.begin(), isas.end());
            terms.emplace_back(tmp);

            //清空数据
            id      =   "";
            name    =   "";
            ns      =   Name_Space::UNKNOWN;
            obs     =   false;
            parts.clear();
            isas.clear();
        }
        if (regex_match(line, id_reg))
        {
            id  = "GO" + line.substr(7);
        }
        if (regex_match(line, name_reg))
        {
            name = line.substr(6);
        }
        if (regex_match(line, ns_reg))
        {
            if (line == "namespace: molecular_function")
            {
                ns = Name_Space::MF;
            }
            
            else if (line == "namespace: cellular_component")
            {
                ns = Name_Space::CC;
            }
            else if (line == "namespace: biological_process")
            {
                ns = Name_Space::BP;
            }
            else
            {
                ns = Name_Space::UNKNOWN;    
            }
            
        }

        if (regex_match(line, obs_reg))
        {
            obs = true;
        }

        if (regex_match(line, part_reg))
        {
            parts.emplace("GO" + line.substr(26,33));
        }
        if (regex_match(line, isa_reg))
        {
            isas.emplace("GO" + line.substr(10,17));
        }

    }

    //加入最后一个节点，
    Term tmp{id,ns,obs,name};

    tmp.get_part_ids().insert(parts.begin(), parts.end());
    tmp.get_isa_ids().insert(isas.begin(), isas.end());
    terms.emplace_back(tmp);


    //移除掉第一个无效节点。
    terms.pop_front();
}

void read_gaf_file(const string file, vector<Annotation> &annos)
{
    ifstream input_file;

    input_file.open(file);

    string line;

    while (getline(input_file, line))
    {
        //获取所有的tab符位置，用于抽取对应的数据
        vector<int> tab_index;
        std::size_t index = 0;
        while ((index =  line.find('\t',index+1)) < (line.size() -1))
        {
            tab_index.push_back(index);
        }

        string go_id;
        string gene_name;
        string evidence_code;
        Name_Space name_space = Name_Space::UNKNOWN;
        set<string> synonym;
        if (tab_index.size() == 16)
        {
            go_id = line.substr(tab_index[3] + 1, 10);
            gene_name = line.substr(tab_index[1] + 1, tab_index[2] - tab_index[1]);
            evidence_code = line.substr(tab_index[5] + 1, tab_index[6] - tab_index[5]);
            string ns = line.substr(tab_index[7] + 1, 1);
            if (ns == "P")
            {
                name_space = Name_Space::BP;
            }
            else if (ns == "F")
            {
                name_space = Name_Space::MF;
            }
            else if (ns == "C")
            {
                name_space = Name_Space::CC;
            }

            string syn = line.substr(tab_index[9] + 1, tab_index[10] - tab_index[9]);
            int index;
            while((index = syn.find("|")) < syn.size())
            {
                synonym.emplace(syn.substr(0,index));

                syn = syn.replace(syn.begin(), syn.begin() + index +1, "");
            }

            Annotation tmp{go_id,gene_name,evidence_code,name_space};
            tmp.get_synonym_gene().insert(synonym.begin(),synonym.end());
            annos.push_back(tmp);
        }
    }

    
}