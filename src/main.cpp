#include "../include/sim.h"

#include <fstream>
#include <iostream>
#include <cassert>
using namespace std;


int main()
{
        // ifstream input;

        // input.open("./result/ids.result");

        // assert(input.is_open());

        // string line;
        // set<string> ids;
        // while (getline(input, line))
        // {
        //         ids.emplace(line);
        // }

        // std::cout << ids.size() << endl;

        // ofstream out("./result/id.result");
        // for (auto line : ids)
        // {
        //         out << line << endl;
        // }
        // out.close();

        // 前三步的计算都能够通过多线程来进行加速
        // 读取参数，选择要进行的计算，输出文件如果已经存在，则不会进行计算，避免覆盖掉已有数据


        //矩阵完备化 MatrixComple
        //输入参数,  类型选择，输出文件，参数列表(字符串，逗号分隔多个参数)，参数解析由完备化的具体类进行解析， 线程数，默认为2
        // 只依赖于net数据

        // 术语相似度计算 TermSimCalc
        // 输入参数， 网络数据文件， 输出文件，线程数，默认为2
        // 依赖于网络数据，本体结构，注释信息
        
        Calculator::TermSim::calculator("./data/net.txt", "./result/term.result", 2);


        // 基因相似度计算   GeneSim
        // 输入参数， 要计算的基因对文件，输出文件，线程数，默认为2
        // 依赖于术语相似度
        // Calculator::LFCValue::gene_pair_generator("./result/pair.result");
        
        // Calculator::GeneSim::calculator("./result/pair.result", "./result/gene.result");

        // lfc计算  Evaluator
        // 输入参数， 基因相似度文件，输出文件，是否输出到控制台
        // 依赖于基因相似度
        
        
        // Calculator::LFCValue::calculator("./result/gene.result", "./result/lfc.result", true);

        return 0;
}