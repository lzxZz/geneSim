#include <iostream>
//#include "../include/shared.h"
#include "../include/data.h"
#include "../include/getter.h"

using namespace std;

int
main(int argc, char **argv){
    Data::Getter getter;
    const set<string>  ss = getter.get_public_ancestor_by_id("GO0015422", "GO0015423");

    cout << ">>>" << ss.size() << "<<<" << endl;
    
    return 0;
}
