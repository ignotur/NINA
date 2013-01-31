#include <iostream>
#include <fstream>
#include "stars.h"

using namespace std;

P_distr::P_distr (ifstream * in) {
    *in>>ain;
    *in>>bin;
}

void P_distr::print_param (ostream * out) {
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// Начальное распределение P и B читается из входных файлов //"<<endl;
    *out<<"// P="<<a()<<", dP="<<b()<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

double P_distr::a () {
    return ain;
}

double P_distr::b () {
    return bin;
}

B_distr::B_distr (ifstream * in) {
    *in>>ain;
    *in>>bin;
}

void B_distr::print_param (ostream * out) {
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// Начальное распределение P и B читается из входных файлов //"<<endl;
    *out<<"// B="<<a()<<", dB="<<b()<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

double B_distr::a () {
    return ain;
}

double B_distr::b () {
    return bin;
}
