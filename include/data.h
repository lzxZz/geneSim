#ifndef __DATA_H
#define __DATA_H

#include <vector>
#include "edge.h"
#include "term.h"
#include "anno.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <deque>


using std::vector;
using std::deque;

using std::ifstream;
using std::istringstream;

//从文件读取net数据，通过vector引用返回数据。
void read_net_file(const string file,vector<Edge> &edges);
void read_obo_file(const string file,deque<Term> &terms);
void read_gaf_file(const string file, vector<Annotation> &annos);
void read_ec_file(const string file, vector<string> &ecs);

vector<string> string_split_by_char(const string &str, const char symbol);




#endif