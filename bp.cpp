#include <iostream>
#include <fstream>
#include "stars.h"

using namespace std;

PDistr::PDistr ()	       {
    ain = 0.3;
    bin = 0.15;
}

void PDistr::Set (ifstream * in) {
    *in>>ain;
    *in>>bin;
}

void PDistr::Set (double avalue, double bvalue)	{
    ain=avalue;
    bin=bvalue;
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

BDistr::BDistr ()	       {
    ain = 12.65;
    bin = 0.55;
}

void BDistr::Set (ifstream * in) {
    *in>>ain;
    *in>>bin;
}

void BDistr::Set (double avalue, double bvalue)	{
    ain=avalue;
    bin=bvalue;
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
