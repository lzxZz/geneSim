#include "../include/data.h"
#include <boost/algorithm/string.hpp>
#include <regex>
#include <fstream>
#include <iostream>
#include <assert.h>

using boost::is_any_of;
using std::ifstream;

set<string> Data::Data::null_set{};

void Data::Data::init_root_count()
{
    if (is_init_root_count)
    {
        return;
    }
    Log::log({__FILE__, __func__, "开始初始化根节点类别下所有的注释基因总数"});
    //初始化注释文件项列表，以获取所有的基因
    init_gaf_list();

    for (Annotation anno : gaf_items)
    {
        switch (anno.get_name_space())
        {
        case Name_Space::BP:
            anno_gene_count_bp++;
            break;
        case Name_Space::CC:
            anno_gene_count_cc++;
            break;
        case Name_Space::MF:
            anno_gene_count_mf++;
            break;
        case Name_Space::UNKNOWN:;
        }
    }
    is_init_root_count = true;
    Log::log({__FILE__, __func__, "root_count 初始化完成"});
}

void Data::Data::init_gaf_list()
{
    if (is_init_gaf_file)
    {
        return;
    }
    Log::log({__FILE__, __func__, "开始初始化注释基因列表"});

    string file = "./data/gene.gaf";

    ifstream input_file;

    input_file.open(file);

    string line;

    while (getline(input_file, line))
    {

        vector<string> infos;
        boost::split(infos, line, boost::is_any_of("\t"));

        //正常的注释文件有17列
        assert(infos.size() == 17);

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
        vector<string> tmp_synonym;
        boost::split(tmp_synonym, infos[10], boost::is_any_of("|"));

        Annotation tmp{go_id, gene_name, evidence_code, name_space};
        tmp.get_synonym_gene().insert(tmp_synonym.begin(), tmp_synonym.end());
        gaf_items.push_back(tmp);
    }
    is_init_gaf_file = true;
    Log::log({__FILE__, __func__, "注释基因列表初始化完成"});
}

void Data::Data::init_obo_list()
{
    if (is_init_obo_file)
    {
        return;
    }
    Log::log({__FILE__, __func__, "开始初始化本体节点列表"});

    string obo_buf = "./data/obo.buf";

    ifstream buf_reader;
    buf_reader.open(obo_buf);

    if (buf_reader.is_open())
    {
        string line;
        while (getline(buf_reader, line))
        {
            vector<string> values;
            boost::split(values, line, boost::is_any_of("|"));
            assert(values.size() == 6);
            string id = values[0];
            string name = values[1];
            Name_Space ns;
            bool obs;
            set<string> iss;
            set<string> parts;

            if (values[2] == "B")
                ns = Name_Space::BP;
            else if (values[2] == "M")
                ns = Name_Space::MF;
            else if (values[2] == "C")
                ns = Name_Space::CC;
            else
                ns = Name_Space::UNKNOWN;

            if (values[3] == "true")
            {
                obs = true;
            }
            else
            {
                obs = false;
            }

            Term term{id, ns, obs, name};
            set<string> tmp_set;
            boost::split(tmp_set, values[4], boost::is_any_of("\t"));
            term.get_isa_ids().insert(tmp_set.begin(), tmp_set.end());
            term.get_isa_ids().erase("");

            boost::split(tmp_set, values[5], boost::is_any_of("\t"), boost::token_compress_on);
            term.get_part_ids().insert(tmp_set.begin(), tmp_set.end());
            term.get_part_ids().erase("");
            onto_items.push_back(term);

            
        }
        buf_reader.close();
    }
    else
    {

        string file = "./data/onto.obo";

        deque<Term> tmp_deque;
        ifstream input_file;
        input_file.open(file);
        if (!input_file.is_open())
        {
            //std::cout << "obo file read errer!\n can't open the file" << std::endl;
        }

        //节点属性参数定义
        string id;
        string name;
        Name_Space ns;
        bool obs;
        set<string> parts;
        set<string> isas;

        regex id_reg("id: GO:\\d{7}");
        regex name_reg("name: .+");
        regex obs_reg("is_obsolete: true");
        regex ns_reg("namespace: .+");
        regex part_reg("relationship: part_of GO:\\d{7}.+");
        regex isa_reg("is_a: GO:\\d{7}.+");

        string line;
        while (getline(input_file, line))
        {
            if (line == "[Term]")
            {
                Term tmp{id, ns, obs, name};

                tmp.get_part_ids().insert(parts.begin(), parts.end());
                tmp.get_isa_ids().insert(isas.begin(), isas.end());
                tmp_deque.emplace_back(tmp);

                //std::cout <<  tmp.debug() << std::endl;

                //清空数据
                id = "";
                name = "";
                ns = Name_Space::UNKNOWN;
                obs = false;
                parts.clear();
                isas.clear();
            }
            if (regex_match(line, id_reg))
            {
                id = "GO" + line.substr(7);
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
                parts.emplace("GO" + line.substr(25, 7));
            }
            if (regex_match(line, isa_reg))
            {
                isas.emplace("GO" + line.substr(9, 7));
            }
        }

        //加入最后一个节点，
        Term tmp{id, ns, obs, name};

        tmp.get_part_ids().insert(parts.begin(), parts.end());
        tmp.get_isa_ids().insert(isas.begin(), isas.end());
        tmp_deque.emplace_back(tmp);

        //移除掉第一个无效节点。
        tmp_deque.pop_front();
        onto_items.insert(onto_items.begin(), tmp_deque.begin(), tmp_deque.end());

        ofstream buf_writer;
        buf_writer.open(obo_buf);
        assert(buf_writer.is_open());

        for (auto term : onto_items)
        {
            buf_writer << term.get_id() << "|";
            buf_writer << term.get_name() << "|";
            buf_writer << static_cast<char>(term.get_name_space()) << "|";
            if (term.is_obsolete())
            {
                buf_writer << "true"
                           << "|";
            }
            else
            {
                buf_writer << "false"
                           << "|";
            }

            for (auto id : term.get_isa_ids())
            {
                buf_writer << id << "\t";
            }
            buf_writer << "|";
            for (auto id : term.get_part_ids())
            {
                buf_writer << id << "\t";
            }
            buf_writer << endl;
            Log::log({__FILE__, __func__, "写入缓存完成"});
        }
    }
    is_init_obo_file = true;
    Log::log({__FILE__, __func__, "初始化本体节点列表完成"});
}

void Data::Data::init_net_map()
{
    if (is_init_net_file)
    {
        return;
    }
    Log::log({__FILE__, __func__, "开始初始化功能网络map数据"});
    string file = "./data/net.txt";
    ifstream input_file;

    input_file.open(file);

    string line;
    while (getline(input_file, line))
    {
        istringstream is(line);
        string g1, g2;
        double weight;
        is >> g1;
        is >> g2;
        is >> weight;
        //edges.push_back(Edge{g1,g2,weight});
        string key = g1 + ":" + g2;
        net_value.emplace(make_pair(key, weight));
    }
    is_init_net_file = true;
    Log::log({__FILE__, __func__, "基因功能网络map数据初始化完成"});
}

void Data::Data::init_ec_list_and_map()
{
    if (is_init_ec_file)
    {
        return;
    }
    Log::log({__FILE__, __func__, "开始初始化生物过程信息数据"});

    string file = "./data/ec.tab";

    ifstream input_file;
    input_file.open(file);
    if (!input_file.is_open())
    {
        return;
    }

    string line;
    set<string> tmp_numbers;
    while (getline(input_file, line))
    {
        vector<string> infos;

        //最后一个选项用于合并多个空白，不开启
        boost::split(infos, line, boost::is_any_of("\t"));

        //正常的ec信息将会转化为5个项
        assert(infos.size() == 5);

        if (infos[2] == "")
        {
            continue;
        }

        tmp_numbers.emplace(infos[2]);

        auto it = ecs_genes.find(infos[2]);
        if (it != ecs_genes.end())
        {
            if (infos[3] != "")
            {
                it->second.emplace(infos[3]);
            }
        }
        else
        {
            if (infos[3] != "")
            {
                ecs_genes.emplace(make_pair(infos[2], set<string>{infos[3]}));
            }
        }
    }

    ec_numbers.insert(ec_numbers.end(), tmp_numbers.begin(), tmp_numbers.end());

    is_init_ec_file = true;
    Log::log({__FILE__, __func__, "生物过程数据初始化完成"});
}

void Data::Data::init_anno_map()
{
    if (is_init_anno_map)
    {
        return;
    }

    Log::log({__FILE__, __func__, "开始初始化id-genes,gene-ids的注释map数据"});

    init_gaf_list();

    assert(is_init_gaf_file);
    assert(gaf_items.size() != 0);

    for (Annotation anno : gaf_items)
    {
        string id = anno.get_go_id();
        string gene_name = anno.get_gene_name();

        //插入id-gene的key-value
        auto id_iter = id_gene_anno.find(id);
        if (id_iter == id_gene_anno.end())
        {
            set<string> tmp_set;
            tmp_set.emplace(gene_name);

            id_gene_anno.emplace(std::make_pair(id, tmp_set));
        }
        else
        {
            id_iter->second.emplace(gene_name);
        }

        auto name_iter = gene_id_anno.find(gene_name);
        if (name_iter == gene_id_anno.end())
        {
            set<string> tmp_set;
            tmp_set.emplace(id);

            gene_id_anno.emplace(make_pair(gene_name, tmp_set));
        }
        else
        {
            name_iter->second.emplace(id);
        }
    }

    is_init_anno_map = true;
    Log::log({__FILE__, __func__, "id-genes,gene-ids注释map数据初始化完成"});
}

void Data::Data::init_descendant()
{
    if (is_init_descendant)
    {
        return;
    }
    Log::log({__FILE__, __func__, "开始初始化子孙节点数据"});
    //从缓存尝试读取数据
    string buf_file = "./data/descendant.buf";

    ifstream input_file;

    input_file.open(buf_file);

    if (input_file.is_open())
    {
        Log::log({__FILE__, __func__, "缓存数据存在，从缓冲文件中加载数据"});
        //缓存数据存在
        string line, key;
        while (getline(input_file, line))
        {
            vector<string> infos;
            boost::split(infos, line, is_any_of(":"));
            assert(infos.size() > 1);
            key = infos[0];

            set<string> values;
            boost::split(values, infos[1], is_any_of("\t"));

            auto it = id_descendant.find(key);

            if (it == id_descendant.end())
            {
                id_descendant.emplace(make_pair(key, values));
            }
            else
            {
                it->second.insert(values.begin(), values.end());
            }
        }
        input_file.close();
    }
    else
    {
        Log::log({__FILE__, __func__, "缓存数据不存在，开始手动计算数据"});
        //缓存数据不存在
        for (auto term : onto_items)
        {
            set<string> descendant;
            string pid = term.get_id();

            for (auto child_term : onto_items)
            {
                auto parent_set = id_ancestor.find(child_term.get_id());

                if (parent_set->second.size() > 0)
                {
                    if (parent_set->second.find(pid) != parent_set->second.end())
                    {
                        descendant.emplace(child_term.get_id());
                    }
                }
            }

            id_descendant.emplace(make_pair(pid, descendant));
        }

        Log::log({__FILE__, __func__, "将子孙节点map数据缓存到硬盘"});
        //写入缓存文件
        ofstream output_file;

        output_file.open(buf_file, ios_base::app); //追加写入

        for (auto iter : id_descendant)
        {
            output_file << iter.first << ":";

            if (iter.second.size() == 0)
            {
                output_file << endl;
                continue;
            }
            for (string id : iter.second)
            {
                output_file << id << "\t";
            }

            output_file << endl;
        }

        output_file.close();
    }

    is_init_descendant = true;
    Log::log({__FILE__, __func__, "id-子孙节点map数据初始化完成"});
}

void Data::Data::init_ancestor()
{
    if (is_init_ancestor)
    {
        return;
    }
    Log::log({__FILE__, __func__, "开始初始化id到祖先节点数据"});

    //前提条件，需要id_term的键值对数据
    if (!is_init_id_term)
    {
        init_id_term();
    }

    assert(is_init_id_term);

    for (auto it : id_term)
    {

        deque<string> add_list;
        add_list.push_back(it.first);

        set<string> ancestor;
        while (add_list.size() > 0)
        {
            string id = add_list[0];
            auto term_iter = id_term.find(id);
            assert(term_iter != id_term.end());
            Term term = term_iter->second;

            ancestor.insert(term.get_isa_ids().begin(), term.get_isa_ids().end());
            ancestor.insert(term.get_part_ids().begin(), term.get_part_ids().end());

            add_list.insert(add_list.end(), term.get_isa_ids().begin(), term.get_isa_ids().end());
            add_list.insert(add_list.end(), term.get_part_ids().begin(), term.get_part_ids().end());

            add_list.pop_front();
        }
        id_ancestor.emplace(make_pair(it.first, ancestor));
    }

    is_init_ancestor = true;
    Log::log({__FILE__, __func__, "id到祖先节点数据初始化完成"});
}

void Data::Data::init_id_term()
{

    if (is_init_id_term)
    {
        return;
    }

    if (!is_init_obo_file)
    {
        init_obo_list();
    }
    Log::log({__FILE__, __func__, "开始初始化id到节点的map数据"});

    for (auto term : onto_items)
    {
        id_term.emplace((make_pair(term.get_id(), term)));
    }

    is_init_id_term = true;
    Log::log({__FILE__, __func__, "id到节点的map数据初始化完成"});
}

//初始化id-set的键值对
void Data::Data::init_child()
{
    if (is_init_child)
    {
        return;
    }
    Log::log({__FILE__, __func__, "开始初始化id到子节点的map数据"});
    //缓存文件
    string buf_file = "./data/child.buf";

    ifstream input_file;

    input_file.open(buf_file);

    if (input_file.is_open())
    { //成功打开缓存文件
        Log::log({__FILE__, __func__, "从缓存中加载数据"});
        string line;
        while (getline(input_file, line))
        {
            string key;
            vector<string> infos;
            boost::split(infos, line, boost::is_any_of(":"));

            assert(infos.size() == 2);

            key = infos[0];
            set<string> values;
            boost::split(values, infos[1], boost::is_any_of("\t"));

            auto it = id_child.find(key);

            if (it == id_child.end())
            {
                id_child.emplace(make_pair(key, values));
            }
            else
            {
                it->second.insert(values.begin(), values.end());
            }
        }
        input_file.close();
    }
    else
    { //没有缓存文件
        Log::log({__FILE__, __func__, "没有缓存文件，手动计算数据"});
        for (auto term : onto_items)
        {
            string pid = term.get_id();
            set<string> childs;
            for (auto child_term : onto_items)
            {
                if (child_term.get_isa_ids().find(pid) != child_term.get_isa_ids().end())
                {
                    childs.emplace(child_term.get_id());
                }

                if (child_term.get_part_ids().find(pid) != child_term.get_part_ids().end())
                {
                    childs.emplace(child_term.get_id());
                }
            }

            id_child.emplace(make_pair(pid, childs));
        }

        Log::log({__FILE__, __func__, "将缓存数据存盘"});

        //写入缓存文件
        ofstream output_file;

        output_file.open(buf_file, ios_base::app);

        for (auto it : id_child)
        {
            output_file << it.first << ":";

            for (auto child : it.second)
            {
                output_file << child << "\t";
            }

            output_file << endl;
        }

        output_file.close();
    }

    is_init_child = true;
    Log::log({__FILE__, __func__, "id到子节点数据存盘"});
}

//key = id1:id2
void Data::Data::init_path_node()
{
    if (is_init_path_node)
    {
        return;
    }

    Log::log({__FILE__, __func__, "开始初始化路径map数据"});

    string buf_file = "./data/path.buf";

    ifstream input_file;

    input_file.open(buf_file);

    if (input_file.is_open())
    { //有缓存文件
        Log::log({__FILE__, __func__, "读取缓存文件"});

        string line;
        while (getline(input_file, line))
        {
            string key;
            vector<string> infos;
            boost::split(infos, line, boost::is_any_of("|"));

            assert(infos.size() == 2);

            set<string> values;

            boost::split(values, infos[1], boost::is_any_of("\t"));
            path_nodes.emplace(make_pair(key, values));
        }
        input_file.close();
    }
    else
    { //无缓存文件，从路径缓存开始计算路径
        Log::log({__FILE__, __func__, "无缓存文件，从path.pair中开始计算路径信息"});
        ifstream pair_file;
        pair_file.open("./dir/path.pair");

        if (!pair_file.is_open())
        {
            return;
        }

        //计算路径
        string line;
        while (getline(pair_file, line))
        {
            string id1, id2;
            vector<string> keys;

            boost::split(keys, line, boost::is_any_of(":"));
            id1 = keys[0];
            id2 = keys[1];
            set<string> nodes;
            nodes = get_path_way_node(id1, id2);
            path_nodes.emplace(make_pair(line, nodes));
        }

        //写入缓存文件
        Log::log({__FILE__, __func__, "开始写入缓存文件"});
        ofstream output_file;
        output_file.open("./data/path.buf");

        for (auto path : path_nodes)
        {
            output_file << path.first << "|";
            for (auto value : path.second)
            {
                output_file << value << "\t";
            }
            output_file << endl;
        }
    }

    is_init_path_node = true;
    Log::log({__FILE__, __func__, "路径Map信息初始化完成"});
}

set<string> Data::Data::get_path_way_node(string child_id, string parent_id)
{
    set<string> nodes;

    if (!is_init_child)
    {
        init_child();
    }
    if (!is_init_descendant)
    {
        init_descendant();
    }

    assert(is_init_child);
    assert(is_init_descendant);

    auto it = id_child.find(parent_id);
    Log::log({child_id, parent_id});
    assert(it != id_child.end());

    deque<string> search_list;
    search_list.insert(search_list.end(), it->second.begin(), it->second.end());

    while (search_list.size() > 0)
    {

        string key = search_list.front();

        //子孙节点中存在child_id则为路径上的节点，否则不是

        auto set_it = id_descendant.find(key);

        if (set_it == id_descendant.end())
        {
            search_list.pop_front();
            continue;
        }
        assert(set_it != id_descendant.end());

        set<string> des_set = set_it->second;

        //子孙节点集合中包括child_id的才是路径上的节点，添加至返回结果
        if (des_set.find(child_id) != des_set.end())
        {
            nodes.emplace(key);

            auto add_it = id_child.find(key);
            if (add_it != id_child.end())
            {
                search_list.insert(search_list.end(), add_it->second.begin(), add_it->second.end());
            }
        }

        search_list.pop_front();
    }

    return nodes;
}

void Data::Data::init_net_matrix()
{
    if (is_init_net_matrix)
    {
        return;
    }

    set<string> gene_names;
    ifstream input;
    input.open("./data/net.txt");
    assert(input.is_open());

    string line;
    while (getline(input, line))
    {
        istringstream is(line);
        string gene;
        is >> gene;
        gene_names.emplace(gene);
        is >> gene;
        gene_names.emplace(gene);
    }

    int index = 0;

    vector<string> genes;
    genes.insert(genes.end(), gene_names.begin(), gene_names.end());


    for (auto gene : genes)
    {
        gene_name_index.emplace(make_pair(gene, index++));
        net_array.push_back(vector<double>{});
    }

    if (!is_init_net_file)
    {
        init_net_map();
    }
    assert(is_init_net_file);

    int count = genes.size();

    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < count; j++)
        {
            string key = genes[i] + ":" + genes[j];
            net_array.at(i).push_back(get_net_value_by_key(key));
        }

    }

    is_init_net_matrix = true;
}

//****************************************************************
//****************************************************************
//**********************公开返回数据********************************
//****************************************************************
//****************************************************************

vector<std::string> Data::Data::get_ec_numbers()
{
    if (!is_init_ec_file)
    {
        init_ec_list_and_map();
    }

    assert(is_init_ec_file);
    assert(ec_numbers.size() != 0);

    return ec_numbers;
}
set<string> Data::Data::get_gene_set_by_ec_number(string ec_number)
{
    if (!is_init_ec_file)
    {
        init_ec_list_and_map();
    }
    assert(is_init_ec_file);
    assert(ecs_genes.size() != 0);

    auto it = ecs_genes.find(ec_number);

    if (it == ecs_genes.end())
    {
        return set<string>{};
    }

    assert(it != ecs_genes.end());
    return it->second;
}

set<string> Data::Data::get_anno_id_set_by_gene(string gene_name)
{
    if (!is_init_anno_map) //减少函数调用
    {
        init_anno_map();
    }

    assert(is_init_anno_map);
    assert(gene_id_anno.size() != 0);

    auto it = gene_id_anno.find(gene_name);

    if (it == gene_id_anno.end())
    {
        return set<string>{};
    }

    assert(it != gene_id_anno.end());
    return it->second;
}

set<string> Data::Data::get_anno_gene_set_by_id(string id)
{
    if (!is_init_anno_map)
    {
        init_anno_map();
    }

    assert(is_init_anno_map);
    assert(gene_id_anno.size() != 0);

    auto it = id_gene_anno.find(id);

    if (it == id_gene_anno.end())
    {
        return set<string>{};
    }

    assert(it != id_gene_anno.end());

    return it->second;
}
set<string> Data::Data::get_descendant_by_id(string id)
{
    if (!is_init_descendant)
    {
        init_descendant();
    }

    assert(is_init_descendant);
    assert(id_descendant.size() != 0);

    auto it = id_descendant.find(id);

    if (it == id_descendant.end())
    {
        return set<string>{};
    }

    assert(it != id_descendant.end());

    return it->second;
}

set<string> Data::Data::get_child_by_id(string id)
{
    if (!is_init_child)
    {
        init_child();
    }

    assert(is_init_child);

    auto it = id_child.find(id);

    if (it == id_child.end())
    {
        return set<string>{};
    }

    assert(it != id_child.end());

    return it->second;
}

set<string> Data::Data::get_path_term_set_by_id(string term_child, string term_parent)
{
    if (!is_init_path_node)
    {
        init_path_node();
    }

    assert(is_init_path_node);

    string key = term_child + ":" + term_parent;

    auto it = path_nodes.find(key);

    if (it == path_nodes.end())
    {
        return set<string>{};
    }

    assert(it != path_nodes.end());

    return it->second;
}

set<string> Data::Data::get_child_anno_gene_set_by_id(string id)
{
    if (!is_init_anno_map)
    {
        init_anno_map();
    }
    assert(is_init_anno_map);

    if (!is_init_descendant)
    {
        init_descendant();
    }
    assert(is_init_descendant);

    auto it = id_descendant.find(id);

    if (it == id_descendant.end())
    {
        return set<string>{};
    }

    assert(it != id_descendant.end());

    set<string> result;
    for (auto child_id : it->second)
    {
        auto des_iter = id_gene_anno.find(child_id);

        if (des_iter != id_gene_anno.end())
        {
            result.insert(des_iter->second.begin(), des_iter->second.end());
        }
    }

    return result;
}
set<string> Data::Data::get_term_node_anno_gene_set_by_id(string id)
{
    if (!is_init_anno_map)
    {
        init_anno_map();
    }
    assert(is_init_anno_map);

    auto it = id_gene_anno.find(id);

    if (it == id_gene_anno.end())
    {
        return set<string>{};
    }

    assert(it != id_gene_anno.end());

    return it->second;
}

set<string> Data::Data::get_public_ancestor_by_id(string term1, string term2)
{
    if (!is_init_ancestor)
    {
        init_ancestor();
    }
    assert(is_init_ancestor);

    set<string> public_ancestor;

    auto it1 = id_ancestor.find(term1);

    if (it1 == id_ancestor.end())
    {
        return null_set;
    }

    auto it2 = id_ancestor.find(term2);

    if (it2 == id_ancestor.end())
    {
        return null_set;
    }

    for (auto id : it1->second)
    {
        if (it2->second.find(id) != it2->second.end())
        {
            public_ancestor.emplace(id);
        }
    }

    return public_ancestor;
}

double Data::Data::get_net_value_by_key(string key)
{

    if (!is_init_net_file)
    {
        init_net_map();
    }
    assert(is_init_net_file);

    auto it = net_value.find(key);

    if (it != net_value.end())
    {
        return it->second;
    }

    return 0;
}


double Data::Data::get_net_value_by_keys(string gene1, string gene2)
{
    if (! is_init_net_matrix)
    {
        init_net_matrix();
    }

    assert(is_init_net_matrix);

    int index1, index2;

    auto it = gene_name_index.find(gene1);
    if (it == gene_name_index.end())
    {
        return 0;
    }
    index1 = it->second;

    it = gene_name_index.find(gene2);
    if (it == gene_name_index.end())
    {
        return 0;
    }
    index2 = it->second;
    

    return net_array.at(index1).at(index2);
}

int Data::Data::get_root_node_anno_gene_count(Name_Space ns)
{
    if (!is_init_root_count)
    {
        init_root_count();
    }

    assert(is_init_root_count);

    switch (ns)
    {
    case Name_Space::BP:
        return anno_gene_count_bp;
    case Name_Space::MF:
        return anno_gene_count_mf;
    case Name_Space::CC:
        return anno_gene_count_cc;
    case Name_Space::UNKNOWN:
        return 0;
    }

    return 0;
}

Term Data::Data::get_node_by_id(string id)
{
    if (!is_init_id_term)
    {
        init_id_term();
    }
    assert(is_init_id_term);

    return id_term.find(id)->second;
}
