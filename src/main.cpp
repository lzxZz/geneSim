#include "../include/sim_term.h"
#include "../include/sim_gene.h"
#include "../include/sim_lfc.h"

#include <fstream>
#include <iostream>
#include <cassert>
using namespace std;
#include "../include/matrix.h"


int main()
{
        // Matrix::Matrix matrix(2,2,10);
        // matrix.print();
        // matrix.set_value(0,0,1000);
        // matrix.print();
        // matrix.multi(0);
        // matrix.print();
        // Matrix::Matrix::getE(5).print();


        // 前三步的计算(矩阵完备化,术语相似度,基因相似度)都能够通过多线程来进行加速
        // 读取参数，选择要进行的计算，输出文件如果已经存在，则不会进行计算，避免覆盖掉已有数据


        //矩阵完备化 MatrixComple
        //输入参数,  类型选择，输出文件，参数列表(字符串，逗号分隔多个参数)，参数解析由完备化的具体类进行解析， 线程数，默认为2
        // 只依赖于net数据

        // 术语相似度计算 TermSimCalc
        // 输入参数， 网络数据文件， 输出文件，线程数，默认为2
        // 依赖于网络数据，本体结构，注释信息
        // Calculator::GeneSim::general("./result/pair.result", "./result/ids.result");
        // Calculator::TermSim::calculator("./data/net.txt", "./result/term.result", 2);

        // 依赖于距离矩阵,基因名称到索引的映射.
        // Calculator::TermSim::calculator_by_matrix("./dir/result.mat","./dir/map.txt", "./result/term_by_matrix.result");

        // 基因相似度计算   GeneSim
        // 输入参数， 要计算的基因对文件，输出文件，线程数，默认为2
        // 依赖于术语相似度
        // Calculator::LFCValue::gene_pair_generator("./result/pair.result");
        
        // Calculator::GeneSim::calculator("./result/pair.result", "./result/gene_matrix.result","./result/term_by_matrix.result");

        // lfc计算  Evaluator
        // 输入参数， 基因相似度文件，输出文件，是否输出到控制台
        // 依赖于基因相似度
        
        
        Calculator::LFCValue::calculator("./result/gene_matrix.result", "./result/lfc_matrix1.result", true);

        return 0;
}