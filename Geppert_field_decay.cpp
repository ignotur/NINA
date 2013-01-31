#include <cmath>
#include "stars.h"

using namespace std;

//-----------------------------------------------------------------
// Период вращения нейтронной звезды. Формулы для периода с учётом 
// убывания магнитного поля получены в основной работе методом
// интегрирования дифференциального уравнения для электромагнитного
// дипольного излучения и токовых потерь

double neutron_star::get_P(double t/*, parametrs_B * param*/) {
double P_res, I, tmp, alpha, tau_ohm;
I = 2./5. * M*M_sol*pow(R,2);
t = t - tau;

alpha   = paramet_B->get_alpha();
tau_ohm = paramet_B->get_tau_ohm();

double tau_mu = 3.*pow(light_velocity,3)*I/2./pow(B,2)/pow(R,6);

	if (t<=1e6)	{
		tmp = alpha*pow(e, -t/tau_ohm);
		P_res = sqrt(P*P - 8*3.2e7*pi*pi*tau_ohm/tau_mu*(log(abs(tmp-alpha-1))/alpha/alpha-(alpha+1)/(tmp*alpha*alpha-alpha*alpha*alpha-alpha*alpha)) + 8*3.2e7*pi*pi*tau_ohm/tau_mu/alpha/alpha*(alpha+1));	}
	else {
		tmp = get_P(tau+1e6);
		P_res = sqrt(tmp*tmp + 16*pi*pi*pow(R, 6)*pow(get_B(t),2)/3./I/pow(light_velocity, 3)*t/lsec);	}

return P_res;
}

//-------------------------------------------------------------------

double neutron_star::get_incl (double t) {
return i_incl;
}

//-------------------------------------------------------------------
// Первая производная периода. Формулы взяты из исходного дифферен-
// циального уравнения

double neutron_star::get_dot_P (double t/*, parametrs_B * param*/) {
double res, I, alpha, tau_ohm;
alpha   = paramet_B->get_alpha();
tau_ohm = paramet_B->get_tau_ohm();
t = t - tau;
I = 2./5. * M*M_sol*pow(R,2);
double tau_mu = 3.*pow(light_velocity,3)*I/2./pow(B,2)/pow(R,6);
/*	if (t>1e6)
	res = pow(e, -2.*1e6/tau_ohm)/get_P(t+tau)/tau_mu/pow(1.+alpha*(1-pow(e, -1e6/tau_ohm)),2);
	else
	res = pow(e, -2.*t/tau_ohm)/get_P(t+tau)/tau_mu/pow(1.+alpha*(1-pow(e, -t/tau_ohm)),2);
*/

	if (t>1e6)
	res = 4*pi*pi/tau_mu*pow(e, -2.*1e6/tau_ohm)/get_P(t+tau)/pow(1.+alpha*(1.-pow(e,-1e6/tau_ohm)), 2);
	else
	res = 4*pi*pi/tau_mu*pow(e, -2*t/tau_ohm)/get_P(t+tau)/pow(1.+alpha*(1.-pow(e,-t/tau_ohm)), 2);



return res;
}
//-------------------------------------------------------------------


//------------------------------------------------------------------
// Падение магнитного поля в соотвествии со статьёй Pons, 2009 
// Формула получена из личной переписки

double neutron_star::get_B (double T/*, parametrs_B * param*/) {
double res_B;
double alpha   = paramet_B->get_alpha();
double tau_ohm = paramet_B->get_tau_ohm();
T = T - tau;
	if (T<1e6)
		res_B = B * pow(e, -T/tau_ohm) / (1.+alpha*(1.-pow(e, -T/tau_ohm)));
	else
		res_B = B * pow(e, -1e6/tau_ohm) / (1.+alpha*(1.-pow(e, -1e6/tau_ohm)));
return res_B;
}
//------------------------------------------------------------------

parametrs_B::parametrs_B (ifstream * in) {
*in>>alpha;
*in>>tau_ohm;
}
/*
parametrs_B::parametrs_B (parametrs_B * param)	{
alpha   = param->get_alpha();
tau_ohm = param->get_tau_ohm();
}
*/
void parametrs_B::print_description (ostream * out) {
*out<<"// Используется модель убывания магнитного поля описанная в //"<<endl;
*out<<"// работе Pons, et al. 2007. Первый миллион лет идёт экспо- //"<<endl;
*out<<"// ненциальное убывание напряжённости магнитного поля, затем//"<<endl;
*out<<"// напряжённость поля не меняется в ходе эволюции.          //"<<endl;
*out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_B::print_parametrs   (ostream * out) {
*out<<"//           Параметры модели убывания поля                 //"<<endl;
*out<<"//----------------------------------------------------------//"<<endl;
*out<<"// alpha   - "<<alpha<<endl;
*out<<"// tau_ohm - "<<tau_ohm<<endl;
*out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_B::print_short       (ostream * out) {
*out<<"alpha   - "<<alpha<<endl;
*out<<"tau_ohm - "<<tau_ohm<<endl;
}

double parametrs_B::get_alpha (void) {
return alpha;
}

double parametrs_B::get_tau_ohm (void) {
return tau_ohm;
}
