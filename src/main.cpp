#include <iostream>
//#include "../include/shared.h"
#include "../include/data.h"
#include "../include/getter.h"



using namespace std;

int
main(int argc, char **argv){
    Data::Getter getter;
    // const set<string>  ss = getter.get_public_ancestor_by_id("GO0015422", "GO0015423");
    // const set<string> pathnodes = getter.get_path_term_set_by_id("GO0000002","GO0006996");
    // cout << ">>>" << ss.size() << "<<<" << endl;
    // cout << getter.get_anno_gene_set_by_id("GO0000165").size() << endl;
    
    // cout << getter.get_child_anno_gene_set_by_id("GO0008150") .size() << endl;
    
    // cout << getter.get_root_node_anno_gene_count(Name_Space::BP) << endl;

    // cout << getter.get_child_by_id("GO0008150").size() << endl;
    // cout << getter.get_descendant_by_id("GO0008150").size() << endl;
    // cout << getter.get_ec_genes_by_number("4.1.1.1").size() << endl;
    // cout << getter.get_ec_numbers().size() << endl;
    // cout << getter.get_net_value_by_key("YML100W:YMR261C") << endl;
    cout << getter.get_term_node_anno_gene_set_by_id("GO0008150").size() << endl;


    return 0;
}
// #include <cstdio>
// void prograss()
// {
//     int rate = 0; 
//     char bar[102]; 
//     char ch[] = "-\\|/"; 
//     bar[0] = '\0';

//     while(rate <= 100)
//     {
//         printf("[%-100s][%d\%][%c]\r", bar, rate, ch[rate%4]);
//         fflush(stdout);
//         bar[rate] = '=';
//         bar[++rate] = '\0';
//         system("sleep 1");
//     }
//     printf("\n");
// }