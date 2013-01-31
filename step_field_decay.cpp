#include <cmath>
#include "stars.h"

using namespace std;

double NeutronStar::get_B (double t) {
    double res;
    t = t - tau;

    if (t < paramet_B->get_time_I()) {
        res = B;
    } else if ((t >= paramet_B->get_time_I()) && (t <= paramet_B->get_time_II())) {
        res = B / paramet_B->get_step();
    } else if (t > paramet_B->get_time_II()) {
        res = 0;
    }

    return res;
}

double NeutronStar::get_incl (double t) {
    return i_incl;
}

double NeutronStar::get_P(double t) {
    double res, I, t_tmp;
    t = t - tau;
    t_tmp = 0;

    I = 2./5. * M*M_sol*pow(R,2);

    if (t > paramet_B->get_time_I())	{
        t_tmp = t;
        t =  paramet_B->get_time_I();
    }

    res = P * sqrt(1. + 16*t*3.2e7*pow(R, 6)*B*B*pi*pi/(3.*pow(light_velocity,3)*I*P*P));

    if (!(t_tmp)) {
        return res;
    }

    if (t_tmp >  paramet_B->get_time_II()) {
        t = paramet_B->get_time_II();
    } else {
        t = t_tmp - paramet_B->get_time_I();
    }

    res = res * sqrt(1. + 16*t*3.2e7*pow(R, 6)*B*B/pow(paramet_B->get_step(), 2)*pi*pi/(3.*pow(light_velocity,3)*I*res*res));

    return res;
}

double NeutronStar::get_dot_P (double t) {
    double res, I;
    I = 2./5. * M*M_sol*pow(R,2);
    //res = pow(get_B(t)/3.2e19, 2)/get_P(t);
    res = 8*pi*pi*pow(R,6)/3./pow(light_velocity, 3)/I/get_P(t)*pow(get_B(t),2);
    return res;
}

parametrs_B::parametrs_B (ifstream * in) {
    *in>>time_I;
    *in>>time_II;
    *in>>step;
}

void parametrs_B::print_description (ostream * out) {
    *out<<"//Учебная модель со ступенчатой функцией затухания магнитно-//"<<endl;
    *out<<"//го поля. Кроме того существует время за которое поле исче-//"<<endl;
    *out<<"//зает.                                                     //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_B::print_parametrs   (ostream * out) {
    *out<<"//      Параметры модели убывания магнитного поля           //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// time I - "<<time_I<<endl;
    *out<<"// time II- "<<time_II<<endl;
    *out<<"// step   - "<<step<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_B::print_short       (ostream * out) {
    *out<<"time I - "<<time_I<<endl;
    *out<<"time II- "<<time_II<<endl;
    *out<<"step   - "<<step<<endl;
}

double parametrs_B::get_time_I (void) {
    return time_I;
}

double parametrs_B::get_time_II (void) {
    return time_II;
}

double parametrs_B::get_step (void) {
    return step;
}

