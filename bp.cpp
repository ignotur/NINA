#include <iostream>
#include <fstream>
#include "stars.h"

using namespace std;

PDistr::PDistr (ifstream * in) {
    *in>>ain;
    *in>>bin;
}

void PDistr::print_param (ostream * out) {
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// Начальное распределение P и B читается из входных файлов //"<<endl;
    *out<<"// P="<<a()<<", dP="<<b()<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

double PDistr::a () {
    return ain;
}

double PDistr::b () {
    return bin;
}

BDistr::BDistr (ifstream * in) {
    *in>>ain;
    *in>>bin;
}

void BDistr::print_param (ostream * out) {
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// Начальное распределение P и B читается из входных файлов //"<<endl;
    *out<<"// B="<<a()<<", dB="<<b()<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

double BDistr::a () {
    return ain;
}

double BDistr::b () {
    return bin;
}
