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

        vector<string> infos = string_split_by_char(line,'\t');
        string go_id;
        string gene_name;
        string evidence_code;
        Name_Space name_space = Name_Space::UNKNOWN;
        set<string> synonym;

        go_id = "GO" + infos[4].substr(3);
        gene_name = infos[2];
        evidence_code = infos[6];

        string ns = infos[8];
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
        vector<string> tmp_synonym = string_split_by_char(infos[10],'|');
       // synonym.emplace(tmp_synonym.begin(), tmp_synonym.end() );

       

        Annotation tmp{go_id,gene_name,evidence_code,name_space};
        tmp.get_synonym_gene().insert(tmp_synonym.begin(),tmp_synonym.end());
        annos.push_back(tmp);
        
    }

    
}

void read_ec_file(const string file, vector<string> &ecs)
{

    ifstream input_file;
    input_file.open(file);
    if (! input_file.is_open())
    {
        return;
    }
    
    string line;
    while (getline(input_file, line))
    {
        vector<string> infos;
        infos = string_split_by_char(line,'\t');
        ecs.emplace_back(infos[2] + "|" + infos[3]);
    }
    
    
}

vector<string> string_split_by_char(const string &str, const char symbol)
{
    vector<string> result;
    
    string::const_iterator start,end;
    start = str.begin();
    end = str.begin();
    
    while(end != str.end())
    {
        if (*end == symbol)
        {
            result.emplace_back(start,end);
            start = end;
            start++;
        }
        end++;
    } 
    return result;
}
