#include <cmath>
#include "stars.h"

using namespace std;

//double int_xi (double, double, double, double, double, double);

//-----------------------------------------------------------------
// Период вращения нейтронной звезды. Формулы для периода с учётом 
// убывания магнитного поля получены в основной работе методом
// интегрирования дифференциального уравнения для электромагнитного
// дипольного излучения и токовых потерь

double neutron_star::get_P(double t/*, parametrs_B * param*/) {
double P_res, I, tmp, tau_Hall, tau_ohm, t_step, I_dp, k_1, k_2, k_3, k_4, kx_1, kx_2, kx_3, kx_4, xi_res, diff_t;
double param, intpart, fracpart, perm;
t = t - tau;

P_res  = P*P;


	if (t <= 3.5e5)												{
		param = t/1.e4;
		fracpart = modf(param, &intpart);
		t_step = 1.e4;
			for (int i=0; i < intpart; i++) 							{
				k_1 = t_step*3.2e7/2./1.e39*pow(get_B(tau+i*t_step), 2.);
				k_2 = t_step*3.2e7/2./1.e39*pow(get_B(tau+(i+0.5)*t_step), 2.);
				k_3 = t_step*3.2e7/2./1.e39*pow(get_B(tau+(i+0.5)*t_step), 2.);
				k_4 = t_step*3.2e7/2./1.e39*pow(get_B(tau+(i+1.0)*t_step), 2.);

				P_res += 0.1666667 * (k_1 + 2*k_2 + 2*k_3 + k_4); 			}
					
		t_step = fracpart*t_step;
		perm = intpart*t_step;

			k_1 = t_step*3.2e7/2./1.e39*pow(get_B(perm+tau), 2.);
			k_2 = t_step*3.2e7/2./1.e39*pow(get_B(perm+tau+0.5*t_step), 2.);
			k_3 = t_step*3.2e7/2./1.e39*pow(get_B(perm+tau+0.5*t_step), 2.);
			k_4 = t_step*3.2e7/2./1.e39*pow(get_B(perm+tau+1.0*t_step), 2.);
      
			P_res += 0.1666667 * (k_1 + 2*k_2 + 2*k_3 + k_4); 			        	}
	else if (t > 3.5e5)											{
		t_step = 1.e4;
				for (int i=0; i < 35; i++) 						{
				k_1 = t_step*3.2e7/2./1.e39*pow(get_B(tau+i*t_step), 2.);
				k_2 = t_step*3.2e7/2./1.e39*pow(get_B(tau+(i+0.5)*t_step), 2.);
				k_3 = t_step*3.2e7/2./1.e39*pow(get_B(tau+(i+0.5)*t_step), 2.);
				k_4 = t_step*3.2e7/2./1.e39*pow(get_B(tau+(i+1.0)*t_step), 2.);
				P_res += 0.1666667 * (k_1 + 2*k_2 + 2*k_3 + k_4); 			}	
		diff_t = t - 3.5e5;
		P_res += 2./1.e39 * pow(get_B(3.7e5 + tau),2.)*diff_t*3.2e7;					}
													

return sqrt(P_res);
}

//-------------------------------------------------------------------


double neutron_star::get_incl(double t){
double res, I_dp;

return I_dp;
}

//-------------------------------------------------------------------
// Первая производная периода. Формулы взяты из исходного дифферен-
// циального уравнения

double neutron_star::get_dot_P (double t/*, parametrs_B * param*/) {
double res, I, tau_Hall, tau_ohm, I_dp;

	res = pow(get_B(t), 2.)/2e39/get_P(t);

return res;
}
//-------------------------------------------------------------------


//------------------------------------------------------------------
// Падение магнитного поля в соотвествии со статьёй Pons, 2009 
// Формула получена из личной переписки

double neutron_star::get_B (double T/*, parametrs_B * param*/) {
double res_B, B_min;
double tau_Hall   = paramet_B->get_tau_Hall();
double tau_ohm = paramet_B->get_tau_ohm();
T = T - tau;

//	if (T > 350e3)
//		res_B=B*pow(pow(0.034*350e3/1e4, 1.17) + 0.84, -1.);
//	else if (T <= 350e3 && T > 50e3)
//		res_B=B*pow(pow(0.034*(T-50e3)/1e4, 1.17) + 0.84, -1.);
//	else if (T < 50e3)
//		res_B = B;


	res_B = B * exp(-T/5.e6);	

return res_B;
}
//------------------------------------------------------------------

parametrs_B::parametrs_B (ifstream * in) {
*in>>tau_Hall;
*in>>tau_ohm;
}
/*
parametrs_B::parametrs_B (parametrs_B * param)	{
alpha   = param->get_alpha();
tau_ohm = param->get_tau_ohm();
}
*/
void parametrs_B::print_description (ostream * out) {
*out<<"//----------------------------------------------------------//"<<endl;
*out<<"// Используется модель убывания магнитного поля с двумя эк- //"<<endl;
*out<<"// споненциальными компонентами,   tau_Hall и вторая tau_ohm//"<<endl;
*out<<"// так же учитывается убывание угла между магнитным полюсом //"<<endl;
*out<<"// и полюсом вращения.                                      //"<<endl;
*out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_B::print_parametrs   (ostream * out) {
*out<<"//           Параметры модели убывания поля                 //"<<endl;
*out<<"//----------------------------------------------------------//"<<endl;
*out<<"// tau_Hall   - "<<tau_Hall<<endl;
*out<<"// tau_ohm - "<<tau_ohm<<endl;
*out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_B::print_short       (ostream * out) {
*out<<"tau_Hall   - "<<tau_Hall<<endl;
*out<<"tau_ohm - "<<tau_ohm<<endl;
}

double parametrs_B::get_tau_Hall (void) {
return tau_Hall;
}

double parametrs_B::get_tau_ohm (void) {
return tau_ohm;
}
