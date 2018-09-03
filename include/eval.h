#ifndef __EVAL_H
#define __EVAL_H

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <unordered_map>
#include <cmath>
#include "../include/calc.h"

using std::set;
using std::vector;
using std::string;
using std::unordered_map;

double get_diff(string gene, set<string> gene1_set, set<string> gene2_set);
bool is_interact(string ec1,string ec2);
void evaluator();

#endif
