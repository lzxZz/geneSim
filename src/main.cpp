#include <iostream>
#include "../include/eval.h"
#include "../include/calc.h"
#include <boost/timer.hpp>
#include <fstream>


using namespace std;

void calc_3w_gene()
{
    int c = 0;
    string line;
    ifstream input;
    input.open("./data/pair.txt");
    while (getline(input, line))
    {
        if (++c > 30){return;}
        istringstream  is(line);
        string g1,g2;
        is >> g1 >> g2;

        boost::timer timer;
        cout << g1 << "\t" << g2 << "\t" <<  gene_sim(g1, g2) << endl;
        cout << timer.elapsed() << endl;


    }
    
}


int
main(int argc, char **argv){
    // evaluator();
    calc_3w_gene();

    // // cout << term_sim("GO0015422", "GO0015423",{}) <<endl;
    // cout << gene_sim("COX9", "QCR8") << endl;

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