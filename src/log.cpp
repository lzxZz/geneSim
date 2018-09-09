#include "../include/log.h"

#include <iostream>

using std::cout;
using std::endl;

void Log::log(initializer_list<string> msg)
{
    for (auto info : msg)
    {
        cout << info << "\t";
    }

    cout << endl;

}