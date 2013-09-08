#include <cmath>
#include "stars.h"


using namespace std;

double MFDConst::get_B (double t, double B, double i_incl, double P) {
    return B;
}

double MFDConst::get_incl (double t, double B, double i_incl, double P) {
    return i_incl;
}

double MFDConst::get_P(double t, double B, double i_incl, double P) {
    double res, I;
    //I = 2./5. * M*M_sol*pow(R,2);
    I=1e45;

    res = P * sqrt(1. + 16*t*3.2e7*pow(1e6, 6)*B*B*pi*pi/(3.*pow(light_velocity,3)*I*P*P));
    //res = sqrt(P+2*B*B*t*3.2e7/1.e39);
    //cout<<R<<endl;
    return res;
}

double MFDConst::get_dot_P (double t, double B, double i_incl, double P) {
    double res, I;
    //I = 2./5. * M*M_sol*pow(R,2);
    //res = pow(B/3.2e19, 2)/get_P(t);
    I=1e45;

    res = 8*pi*pi*pow(1e6,6)/3./pow(light_velocity, 3)/I/get_P(t, B, i_incl, P)*B*B;
    //res = 1*B*B/1.e39/get_P(t);

    return res;
}

MFDConst::MFDConst (vector <double> * values) {
}

void MFDConst::print_description (ostream * out) {
    *out<<"#//           Модель без убывания магнитного поля            //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

void MFDConst::print_parameters   (ostream * out) {
    *out<<"#//   Так как поле не убывает параметров у модели тоже нет   //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

