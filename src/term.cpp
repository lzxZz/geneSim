#include "../include/sim.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

using std::ofstream;



double Calculator::TermSim::get_term_sim_by_ids(string term1, string term2, initializer_list<string> ignore_genes)
{
    ofstream out;
    out.open("./result/idgenes.result",std::ios_base::app);
    ofstream out1;
    out1.open("./result/ids.result",std::ios_base::app);

    string g1 = *ignore_genes.begin();
    string g2 = *(ignore_genes.begin()+1);
    if (term1 < term2)
    {
        out << term1 << "\t" << term2;
        out1 << term1 << "\t" << term2 << endl;
        if (g1 < g2)
        {
            out << "\t" << g1 << ":" << g2 << endl;
        }
        else
        {
            out << "\t" << g2 << ":" << g1 << endl;
        }

        cout << term1 << "\t" << term2 << endl;
    }
    else
    {

        out << term2 << "\t" << term1;
        out1 << term2 << "\t" << term1 << endl;
        if (g1 < g2)
        {
            out << "\t" << g1 << ":" << g2 << endl;
        }
        else
        {
            out << "\t" << g2 << ":" << g1 << endl;
        }
        cout << term2 << "\t" << term1 << endl;
    }
    return 0;
}