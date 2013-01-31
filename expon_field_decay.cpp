#include <cmath>
#include "stars.h"

using namespace std;

double neutron_star::get_B (double t) {
    //double res, A_1 = 0.85766, A_2 = 0.14234;
    double res, A_1 = 0.65, A_2 = 0.35;

    t = t - tau;

    res = B*A_1*(exp(-t/paramet_B->get_tau_ohm()) + A_2/A_1);

    return res;
}

double neutron_star::get_incl (double t) {
    return i_incl;
}

double neutron_star::get_P(double t) {
    //double res, I, t_tmp, A_1 = 0.85766, A_2 = 0.14234;
    double res, I, t_tmp, A_1 = 0.65, A_2 = 0.35;

    t = t - tau;

    I = 2./5. * M*M_sol*pow(R,2);

    res = sqrt(P*P+16.*pi*pi*pow(R,6)*B*B*A_1*A_1/3./pow(light_velocity, 3)/I*(3.156e7*paramet_B->get_tau_ohm()/2.*(1-exp(-2.*t/paramet_B->get_tau_ohm()))+2*3.156e7*paramet_B->get_tau_ohm()*A_2/A_1*(1-exp(-t/paramet_B->get_tau_ohm())) + 3.156e7*pow(A_2/A_1, 2)*t));

    //double dt = t/50.;
    //res = P*P;

    //	for (int i=0; i < 50; i++)
    //		res += pow(get_B(tau+dt*i)/3.2e19, 2.)*dt*2e7;

    //res = sqrt(res);

    return res;
}

double neutron_star::get_dot_P (double t) {
    double res, I;
    I = 2./5. * M*M_sol*pow(R,2);
    //res = 8*pi*pi*pow(R,6)/3./pow(light_velocity, 3)/I/get_P(t)*pow(get_B(t),2);

    res = pow(get_B(t)/3.2e19, 2.)/get_P(t);

    return res;
}

parametrs_B::parametrs_B (ifstream * in) {
    *in>>tau_ohm;
}

void parametrs_B::print_description (ostream * out) {
    *out<<"//Практическая модель представляющая из себя сумму экспонен-//"<<endl;
    *out<<"//и постоянного значения. Обладает гладкостью, в отличие от //"<<endl;
    *out<<"// второй и третей модели. Подбирается только время         //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_B::print_parametrs   (ostream * out) {
    *out<<"//      Параметры модели убывания магнитного поля           //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// tau_ohm - "<<tau_ohm<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

//void parametrs_B::print_short   (ostream * out) {
//*out<<"tau_ohm - "<<tau_ohm<<endl;
//}

void parametrs_B::print_short       (ostream * out) {
    *out<<"tau_ohm - "<<tau_ohm<<endl;
}


double parametrs_B::get_tau_ohm (void) {
    return tau_ohm;
}
