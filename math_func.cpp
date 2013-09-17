#include <cmath>
#include "stars.h"
#include <cstdlib>
#include <iostream>

using namespace std;

double norm_distr ();

double RDFFaucher::rho (double r) {
    double const a = 1.64, b = 4.01, R_1 = 0.55, R_sol = 8.5;
    double const A = 0.0545821764459; // константа получена численным интегрированием

    double res;

    res = A * pow(((r+R_1)/(R_sol + R_1)), a) * pow(e, (-b*((r-R_sol)/(R_sol+R_1))));

    return res;
}

void  RDFFaucher::print_description (ostream * out) {
*out<<"#The radial distribution function of stars in the Galaxy is as in"<<endl;
*out<<"#article by Faucher-Giguere, Kaspi (2006)"<<endl;
*out<<"#----------------------------------------------------------------"<<endl;
}

//--------------------------------------------//
// Другие варинты радиального распределения

// Распределение из работы van der Kruit, 1987
// на основе распределения наблюдаемой поверхностной
// яркости в полосе J Sc галактик

double RDFKruit::rho (double r) {
    double res, R_exp = 4.5, a_R = 1.0683;
    res = a_R*r/pow(R_exp, 2)*exp(-r/R_exp);
    return res;
}

void  RDFKruit::print_description (ostream * out) {
*out<<"#The radial distribution function of stars in the Galaxy is as in"<<endl;
*out<<"#article by van der Kruit (1987). Based on observed distribution"<<endl;
*out<<"#of brightness in J band for Sc galaxies."<<endl;
*out<<"#----------------------------------------------------------------"<<endl;
}
// B00 модель построенная на основе данных
// о распределении яркости в далёком инфракрасном диапазоне
// и милиметрового излучения

double RDFB0::rho  (double r) {
    double res, r_exp = 1.78, sigma = 2.38, r_centr = 4.7, weight = 3.8781977;

    if (r <= r_centr) {
        res = exp(-pow((r-r_centr),2)/pow(sigma,2)) / weight;
    }

    if (r >  r_centr) {
        res = exp((r_centr - r) / r_exp) / weight;
    }

    return res;
}
void  RDFB0::print_description (ostream * out) {
*out<<"#The radial distribution function of stars in the Galaxy is"<<endl;
*out<<"#based on information about distribution of brightness in "<<endl;
*out<<"#in far IR and milimeter radiation."<<endl;
*out<<"#----------------------------------------------------------------"<<endl;
}
// Модель, построенная по распределению поверхностной плотности
// остатков взрывов сверхновых
double RDFSN::rho (double r) {
    double res, alpha = 2, beta = 3.53, R_0 = 8.5;
    res = pow(r/R_0, alpha)*exp(-beta*(r-R_0)/R_0);
    return res;
}

void  RDFSN::print_description (ostream * out) {
*out<<"#The radial distribution function of stars in the Galaxy is"<<endl;
*out<<"#based on distribution of volume density of SN remnants "<<endl;
*out<<"#----------------------------------------------------------------"<<endl;
}
// Модель построенная на наблюдении пульсаров

double RDFPuls::rho (double r) {
    double res, R_peak = 7.04, sigma = 1.83;
    res = 1/sqrt(2*pi) / sigma * exp(-pow(r-R_peak, 2)/2/pow(sigma,2));
    return res;
}

void  RDFPuls::print_description (ostream * out) {
*out<<"#The radial distribution function of stars in the Galaxy is"<<endl;
*out<<"#based on distribution of pulsars (old) "<<endl;
*out<<"#----------------------------------------------------------------"<<endl;
}

GDGauss::GDGauss (vector <double> * val) {

if (val->size() != 4)	{
	print_error_parameters_not_enough ();
	exit(3);	
}
else
	values = val;
}

double GDGauss::generate_next () {
double res;
	res = values->at(1) + values->at(3) * norm_distr();	
return res;
}

GDMGauss::GDMGauss (vector <double> * val) {

if (val->size() < 5)	{
	print_error_parameters_not_enough ();
	exit(4);	
}
else
	values = val;
}

void GDGauss::print_param (ostream * out){
*out<<"#Parameters of gaussian distribution: center - "<<values->at(1)<<", standard deviation - "<<values->at(3)<<endl;
}

void GDMGauss::print_param (ostream * out){
*out<<"#Parameters of multi-gaussian distribution:"<<endl;
	for (int i=0; i < values->size()/6; i++)
		 *out<<"#weight - "<<values->at(6*i+1)<<", center - "<<values->at(6*i+3)<<", standard deviation - "<<values->at(6*i+5)<<endl;
}

double GDMGauss::generate_next () {
double res;
double chance_1, sum=0;

        chance_1 = rand () / rand_high_board;
	for (int i=0; i < values->size()/6; i++)	{
		sum += values->at(6*i+1);
			if (chance_1 < sum)		{
				res = values->at(6*i+3) + values->at(6*i+5) * norm_distr();	
				return res;
			}
	}
		
return res;
}

double expon_vel(double y) {
    double res, v_l = 180;
    res = -v_l*log(2*v_l*y);
    return res;
}

//-------------------------------------------//


// Генерация нормально распределённой случайной величины
// преобразованием Бокса-Мюлера

double norm_distr () {

    double chance_1, chance_2, s, dx, dy;
    bool   is_position_set = false;

    do {
        chance_1 = rand () / rand_high_board;
        chance_2 = rand () / rand_high_board;
        chance_1 = 2*(chance_1 - 0.5);
        chance_2 = 2*(chance_2 - 0.5);
        s = chance_1*chance_1 + chance_2*chance_2;

        if (s != 0 && s<=1) {
            is_position_set = true;
        }
    } while (!(is_position_set));

    dx = chance_1 * sqrt(-2*log(s)/s);
    //dy = chance_2 * sqrt(-2*log(s)/s);

    return dx;
}

//-------------------------------------------------------------------//
// Функции - частные производные потенциала Галактики
// Потенциал взят из работы Kuijken & Gilmore 1989 по сути
double M_dh = 1.45e+11*M_sol, M_b = 9.3e+9*M_sol, M_n = 1.e+10*M_sol;
double beta_1 = 0.4, beta_2 = 0.5, beta_3 = 0.1;
double h_1 = 0.325, h_2 = 0.090, h_3 = 0.125;
double a_G = 2.4;
double b_dh = 5.5, b_b = 0.25, b_n = 1.5;

double dphi_dx (double x, double y, double z) {
    double res;

    res = (M_dh*x*G)/pow(pow(a_G+beta_3*sqrt(z*z+h_3*h_3)+beta_2*sqrt(z*z+h_2*h_2)+beta_1*sqrt(z*z+h_1*h_1),2)+y*y+x*x+b_dh*b_dh, 3./2.) + (M_b*x*G)/pow(y*y+x*x+b_b*b_b, 3./2.)+ (M_n*x*G)/pow(y*y+x*x+b_n*b_n, 3./2.);
    return res;
}

double dphi_dy (double x, double y, double z) {
    double res;

    res = (M_dh*y*G)/pow(pow(a_G+beta_3*sqrt(z*z+h_3*h_3)+beta_2*sqrt(z*z+h_2*h_2)+beta_1*sqrt(z*z+h_1*h_1),2)+y*y+x*x+b_dh*b_dh, 3./2.) + (M_b*y*G)/pow(y*y+x*x+b_b*b_b, 3./2.)+ (M_n*y*G)/pow(y*y+x*x+b_n*b_n, 3./2.);
    return res;
}

double dphi_dz (double x, double y, double z) {
    double res;

    res = (M_dh*((beta_3*z)/sqrt(z*z+h_3*h_3)+(beta_2*z)/sqrt(z*z+h_2*h_2)+(beta_1*z)/sqrt(z*z+h_1*h_1))*(a_G+beta_3*sqrt(z*z+h_3*h_3)+beta_2*sqrt(z*z+h_2*h_2)+beta_1*sqrt(z*z+h_1*h_1))*G)/pow(pow(a_G+beta_3*sqrt(z*z+h_3*h_3)+beta_2*sqrt(z*z+h_2*h_2)+beta_1*sqrt(z*z+h_1*h_1),2)+y*y+x*x+b_dh*b_dh, 3./2.);
    return res;
}

double dphi_dr (double x, double y, double z) {
    double res;
    double r = sqrt(x*x+y*y);

    res = (M_dh*r*G)/pow(pow(a_G+beta_3*sqrt(z*z+h_3*h_3)+beta_2*sqrt(z*z+h_2*h_2)+beta_1*sqrt(z*z+h_1*h_1),2)+r*r+b_dh*b_dh, 3./2.) + (M_b*r*G)/pow(r*r+b_b*b_b, 3./2.) + (M_n*r*G)/pow(r*r+b_n*b_n, 3./2.);
    return res;
}

double phi (double x, double y, double z) {
    double res;
    double r = sqrt(x*x+y*y);
    res = -(G*M_dh) / (pow(a_G + beta_3*sqrt(z*z+h_3*h_3)+beta_2*sqrt(z*z+h_2*h_2)+beta_1*sqrt(z*z+h_1*h_1), 2)+b_dh*b_dh+r*r) - (M_b*G)/sqrt(b_b*b_b+r*r) - (M_n*G)/sqrt(b_n*b_n+r*r);
    return res;
}
//----------------------------------------------------------------------//

/*
//----------------------------------------------------------------------//
// Новый гравитационный потенциал нашей Галактики. Взят из работы
// Flynn, Sommer-Larsen & Christensen
double r_0 = 8.5, V_H = 220*1e5*lcm/lsec;
double r_C1 = 2.7, r_C2 = 0.42, b = 0.3;
double M_C1 = 3.e9*M_sol, M_C2 = 1.6e10*M_sol, M_D1 = 6.6e10*M_sol, M_D2 = -2.9e10*M_sol, M_D3 = 3.3e9*M_sol;
double a1 = 5.81, a2 = 17.43, a3 = 34.86;

double dphi_dR (double x, double y, double z)	{
double res;
double r = sqrt(x*x+y*y+z*z);
res = r*V_H*V_H/(pow(r_0, 2) + pow(r, 2)) + r*G*M_C2/pow(pow(r_C2,2)+pow(r,2), 1.5)+ r*G*M_C1/pow(pow(r_C1,2)+pow(r,2), 1.5);

return res;
}

double dphi_dr (double x, double y, double z)	{
double res;
double R = sqrt(x*x+y*y);
res = G*R*M_D3/pow(R*R+pow(sqrt(z*z+b*b)+a3, 2), 1.5) + G*R*M_D2/pow(R*R+pow(sqrt(z*z+b*b)+a2, 2), 1.5) + G*R*M_D1/pow(R*R+pow(sqrt(z*z+b*b)+a1, 2), 1.5) + dphi_dR(x,y,z)*2*R;

return res;
}

double dphi_dx (double x, double y, double z)	{
double res;
double R = sqrt(x*x+y*y);
double r = sqrt(x*x+y*y+z*z);

res = dphi_dr(x,y,z)*x/sqrt(x*x+y*y) + dphi_dR(x,y,z)*x/sqrt(z*z+y*y+x*x);

return res;
}

double dphi_dy (double x, double y, double z)	{
double res;
double R = sqrt(x*x+y*y);
double r = sqrt(x*x+y*y+z*z);

res = dphi_dr(x,y,z)*y/sqrt(x*x+y*y) + dphi_dR(x,y,z)*y/sqrt(z*z+y*y+x*x);

return res;
}

double dphi_dz (double x, double y, double z)	{
double res;
double R = sqrt(x*x+y*y);
double r = sqrt(x*x+y*y+z*z);

res = dphi_dR(x,y,z)*z/sqrt(z*z+y*y+x*x);

return res;
}

//-------------------------------------------------------------------------//
*/
//--------------------------------------------------------//
// Функция задающия конкретный вид правой части
// дифференциального уравнения
// Решается система уравнение вида \dot s = f(s)
// где s = {v_x, v_y, v_z, x, y, z}
// а f(s)= {-\nabla phi_g, v_x, v_y, v_z}
// начальные условия s (tau) = {x,y, z, v_x, v_y, v_z}

void diff_equi (int n, double * input) {
    double result [6];
    result[3] = - dphi_dx(input[0], input[1], input[2]);
    result[4] = - dphi_dy(input[0], input[1], input[2]);
    result[5] = - dphi_dz(input[0], input[1], input[2]);
    result[0] = input [3];
    result[1] = input [4];
    result[2] = input [5];

    for (int i = 0; i < n; i++) {
        input[i] = result[i];
    }
}

// Функция плотности начальной функции масс
// массивных звёзд. Взята из работы P. Kroupa, 2007
// нормирована из условия что M in [6,10]M_sol

double rho_m (double m) {
    return 24.7672624 * pow(m, -2.35);
}

// Функция плотности распределения пульсаров по скоростям Хартмана
// по ссылке научного руководителя на статью Hartman, 1997

double rho_hartman (double u) {
    double res;
    res = 4*pi/pow(1+pow(u,2), 2);
    return res;
}

