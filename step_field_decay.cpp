#include <cmath>
#include <gsl/gsl_vector.h>
#include "stars.h"

using namespace std;

double MFDStep::get_B (double t, double B, double i_incl, double P) {
    double res;

    if (t < t_1) {
        res = B;
    } else if ((t >= t_1) && (t <= t_2)) {
        res = B / step;
    } else if (t > t_2) {
        res = 0;
    }

    return res;
}
double MFDStep::get_incl (double t, double B, double i_incl, double P) {
    return i_incl;
}
double MFDStep::get_P(double t, double B, double i_incl, double P) {
    double res, I, R = 1e6, t_tmp;
    t_tmp = 0;

      I = 1e45;

    if (t > t_1)	{
        t_tmp = t;
        t =  t_1;
    }

    res = P * sqrt(1. + 16*t*3.2e7*pow(R, 6)*B*B*pi*pi/(3.*pow(light_velocity,3)*I*P*P));

    if (!(t_tmp)) {
        return res;
    }

    if (t_tmp > t_2) {
        t = t_2;
    } else {
        t = t_tmp - t_1;
    }

    res = res * sqrt(1. + 16*t*3.2e7*pow(R, 6)*B*B/pow(step, 2)*pi*pi/(3.*pow(light_velocity,3)*I*res*res));

    return res;
}

double MFDStep::get_dot_P (double t, double B, double i_incl, double P) {
    double res, I, R;
    I = 1e45; R = 1e6;
    res = 8*pi*pi*pow(R,6)/3./pow(light_velocity, 3)/I/get_P(t, B, i_incl, P)*pow(get_B(t, B, i_incl, P),2);
    return res;
}

MFDStep::MFDStep (vector <double> * values) {
	t_1  = values->at(1);
	t_2  = values->at(3);
	step = values->at(5);

}

void MFDStep::print_description (ostream * out) {
    *out<<"#//Model with step-like magnetic field decay. It has two time//"<<endl;
    *out<<"#//one for first decay and the second for vanish of the field//"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

void MFDStep::print_parameters   (ostream * out) {
    *out<<"#//      Parameters of magnetic field decay model            //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
    *out<<"#// time I - "<<t_1<<endl;
    *out<<"#// time II- "<<t_2<<endl;
    *out<<"#// step   - "<<step<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

