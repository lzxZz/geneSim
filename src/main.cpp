#include <iostream>
#include "../include/data.h"
#include "../include/term.h"
#include "../include/anno.h"
#include <dirent.h>

using namespace std;

int
main(int argc, char **argv){
    vector<Edge> edges;
    deque<Term>   terms;
    vector<Annotation> annos;
    vector<string> ecs;
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
    

    
    return 0;
}