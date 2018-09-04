#include <iostream>
#include "../include/shared.h"

using namespace std;

extern unordered_map<string, double> net_value;
int
main(int argc, char **argv){
    // vector<Edge> edges;
    // deque<Term>   terms;
    // vector<Annotation> annos;
    // vector<string> ecs;
    // read_net_file("./data/net.txt",edges);
    // read_obo_file("./data/onto.obo",terms);
    // read_gaf_file("./data/gene.gaf",annos);
    // read_ec_file("./data/ec.tab", ecs);
    // cout << edges.size() << endl;
    // cout << terms.size() << endl;
    
    // for (int i = 0; i < 10; i++)
    // {
    //     // cout << terms.at(i).debug() << endl;
    //     cout << ecs.at(i) << endl;
    // }
    init_net_value();
    
    init_ec_tab();
    init_term();

    cout << id_term.size() << endl;
    
    
    return 0;
}