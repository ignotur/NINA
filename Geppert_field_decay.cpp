#include <cmath>
#include "stars.h"

using namespace std;

//-----------------------------------------------------------------
// Период вращения нейтронной звезды. Формулы для периода с учётом
// убывания магнитного поля получены в основной работе методом
// интегрирования дифференциального уравнения для электромагнитного
// дипольного излучения и токовых потерь
double MFDPons::get_P (double t, double B, double i_incl, double P) {
    double P_res, I, tmp;
//    I = 2./5. * M*M_sol*pow(R,2);
//    t = t - tau;
    I = 1e45;
//    alpha   = paramet_B->get_alpha();
//    tau_ohm = paramet_B->get_tau_ohm();

    double tau_mu = 3.*pow(light_velocity,3)*I/2./pow(B,2)/pow(1e6,6);

    if (t<=1e6)	{
        tmp = alpha*pow(e, -t/tau_ohm);
        P_res = sqrt(P*P - 8*3.2e7*pi*pi*tau_ohm/tau_mu*(log(abs(tmp-alpha-1))/alpha/alpha-(alpha+1)/(tmp*alpha*alpha-alpha*alpha*alpha-alpha*alpha)) + 8*3.2e7*pi*pi*tau_ohm/tau_mu/alpha/alpha*(alpha+1));
    } else {
        tmp = get_P(1e6, B, i_incl, P);
        P_res = sqrt(tmp*tmp + 16*pi*pi*pow(1e6, 6)*pow(get_B(t, B, i_incl, P),2)/3./I/pow(light_velocity, 3)*t/lsec);
    }

    return P_res;
}

//-------------------------------------------------------------------
double MFDPons::get_incl (double t, double B, double i_incl, double P) {
    return i_incl;
}

//-------------------------------------------------------------------
// Первая производная периода. Формулы взяты из исходного дифферен-
// циального уравнения
double MFDPons::get_dot_P (double t, double B, double i_incl, double P) {
    double res, I;
//    alpha   = paramet_B->get_alpha();
//    tau_ohm = paramet_B->get_tau_ohm();
//    t = t - tau;
//    I = 2./5. * M*M_sol*pow(R,2);
    I=1e45;
    double tau_mu = 3.*pow(light_velocity,3)*I/2./pow(B,2)/pow(1e6,6);

    if (t>1e6) {
        res = 4*pi*pi/tau_mu*pow(e, -2.*1e6/tau_ohm)/get_P(t, B, i_incl, P)/pow(1.+alpha*(1.-pow(e,-1e6/tau_ohm)), 2);
    } else {
        res = 4*pi*pi/tau_mu*pow(e, -2*t/tau_ohm)/get_P(t, B, i_incl, P)/pow(1.+alpha*(1.-pow(e,-t/tau_ohm)), 2);
    }

    return res;
}
//-------------------------------------------------------------------


//------------------------------------------------------------------
// Падение магнитного поля в соотвествии со статьёй Pons, 2009
// Формула получена из личной переписки

double MFDPons::get_B (double T, double B, double i_incl, double P) {
    double res_B;
//    double alpha   = paramet_B->get_alpha();
//    double tau_ohm = paramet_B->get_tau_ohm();
//    T = T - tau;

    if (T<1e6) {
        res_B = B * pow(e, -T/tau_ohm) / (1.+alpha*(1.-pow(e, -T/tau_ohm)));
    } else {
        res_B = B * pow(e, -1e6/tau_ohm) / (1.+alpha*(1.-pow(e, -1e6/tau_ohm)));
    }

    return res_B;
}
//------------------------------------------------------------------

void MFDPons::print_description (ostream * out) {
    *out<<"#// Model of magnetic field decay is the same as in article  //"<<endl;
    *out<<"#// by Pons et al. (2007). Before 1e6 years field decays and //"<<endl;
    *out<<"#// after it is constant.                                    //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

void MFDPons::print_parameters   (ostream * out) {
    *out<<"//           Parameters of magnetic field decay             //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// alpha   - "<<alpha<<endl;
    *out<<"// tau_ohm - "<<tau_ohm<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

MFDPons::MFDPons (vector <double> * values) {
alpha   = values->at(1);
tau_ohm = values->at(3);
}

