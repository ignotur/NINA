#include <cmath>
#include "stars.h"

using namespace std;
double MFDExpon::get_B (double t, double B, double i_incl, double P) {
    double res, A_1 = 0.65, A_2 = 0.35;

//    t = t - tau;

    res = B*(exp(-t/tau_ohm);

    return res;
}

double MFDExpon::get_incl (double t, double B, double i_incl, double P) {
    return i_incl;
}

double MFDExpon::get_P(double t, double B, double i_incl, double P) {
    //double res, I, t_tmp, A_1 = 0.85766, A_2 = 0.14234;
    double res, I, t_tmp, A_1 = 0.65, A_2 = 0.35;

    res = sqrt(P*P + 1.6e-39 * B*B * 3.156e7 * tau_ohm * (1.0 - exp(-2.0 * t / tau_ohm)));


//      I = 1e45;

//    res = sqrt(P*P+16.*pi*pi*pow(1e6,6)*B*B*A_1*A_1/3./pow(light_velocity, 3)/I*(3.156e7*tau_ohm/2.*(1-exp(-2.*t/tau_ohm))+2*3.156e7*tau_ohm*A_2/A_1*(1-exp(-t/tau_ohm)) + 3.156e7*pow(A_2/A_1, 2)*t));

    return res;
}

double MFDExpon::get_dot_P (double t, double B, double i_incl, double P) {
    double res, I;
    I = 1e45;
    
    res = 1.6e-39 * pow(get_B(t, B, i_incl, P), 2.)/get_P(t, B, i_incl, P);

    return res;
}

MFDExpon::MFDExpon (vector <double> * values) {
tau_ohm = values->at(1);
}

void MFDExpon::print_description (ostream * out) {
    *out<<"#//This model is the sum of exponential and constant value   //"<<endl;
    *out<<"#//Magnetic field decay law is smooth.                       //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

void MFDExpon::print_parameters   (ostream * out) {
    *out<<"#//      Parameters of magnetic field decay                  //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
    *out<<"#// tau_ohm - "<<tau_ohm<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}


