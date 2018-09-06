#include <iostream>
//#include "../include/shared.h"
#include "../include/data.h"
#include "../include/getter.h"

using namespace std;

int
main(int argc, char **argv){
    Data::Getter getter;
    getter.get_public_ancestor_by_id("","");
    
    return 0;
}
