#include <cstdlib>
#include <time.h>
#include <cmath>
#include <iostream>
#include "stars.h"
#include <algorithm>
#include <cstring>

using namespace std;

//-------------------------------------------------------------------------------------//
// neutron_star --- a class which create neutron star on the place where massive star
// exploded. This class also allows to trace its following evolution. The constructor
// should get an exemplar of the class massive star. The position of the star, 
// magnetic flux and mass of the core are used.
// neutron_star --- класс, который на месте массивной звезды создаёт нейтронную звезду,
// а также позволяет прослеживать её дальнейшую эволюцию. Конструктору класса передаётся
// звезда предшественник этой нейтронной звезды. Извлекаеться положение в пространстве
// Галактике, магнитный поток, масса ядра прародителя
// Масса по статье Hurley, 2000
//
//  Author: Igoshev Andrei
//  Adviser: Alexander F. Kholtygin
//  e-mail: ignotur@gmail.com
//-------------------------------------------------------------------------------------//

//------------------------------------------------------------------//
// Декларация функций
double rho_hartman (double);
double norm_distr  (void);
double expon_vel   (double);
void   Runge_Kutta (int, double, double *, void (*f)(int, double *));
void   diff_equi   (int, double *);
//------------------------------------------------------------------//

NeutronStar::NeutronStar (double T, OBStar  * proteg, MFD * mfd, LM * lm, GD *p_init, GD *b_init) {
    double v, prover, sigma, phi, psi, f, chance_1, chance_2, chance_3;
    bool is_position_set = false;

    M = 1.17 + 0.09 * proteg->get_M_c_SN();

    // Determination of NS radius -------------------
    // according to Latimer et al. (2000)
    R = 3.04*G_cgs*M*M_sol/pow(light_velocity, 2);
    R = R / sqrt(1-2*G_cgs*M*M_sol/(R*pow(light_velocity,2)));
    //-----------------------------------------------------------

    // Generate initial period of pulsar
    do {

    P = p_init->generate_next();
    } while (P < 0);

    // The initial position of the NS is the same as massive star had
    x = proteg->get_position_x();
    y = proteg->get_position_y();
    z = proteg->get_position_z();

    //-------------------------------------------------------------
    // Here the natal kick is imparted 
    // The values and prescription for the natal kick is from the article by Verbunt, Igoshev & Cator (2017)
    double is_velocity_set;

    chance_1 = rand () / rand_high_board;
    
    if (chance_1 < 0.42) {

	v_x = 75.0*1e5*lcm/lsec*norm_distr();	
 	v_y = 75.0*1e5*lcm/lsec*norm_distr();	
	v_z = 75.0*1e5*lcm/lsec*norm_distr();	
    }
    else {
        v_x = 316.0*1e5*lcm/lsec*norm_distr();   
        v_y = 316.0*1e5*lcm/lsec*norm_distr();
        v_z = 316.0*1e5*lcm/lsec*norm_distr();
    }

    /*	
    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_x = 130.0*1e5*lcm/lsec*norm_distr();
    } else {
        v_x = -130.0*1e5*lcm/lsec*norm_distr();
    }

    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_y = 130.0*1e5*lcm/lsec*norm_distr();
    } else {
        v_y = -130.0*1e5*lcm/lsec*norm_distr();
    }

    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_z = 60.0*1e5*lcm/lsec*norm_distr();
    } else {
        v_z = -60.0*1e5*lcm/lsec*norm_distr();
    }
    */ 
    //v_x =  proteg->get_velocity_x();
    //v_y =  proteg->get_velocity_y();
    //v_z =  proteg->get_velocity_z();



    B=b_init->generate_next();
    B=pow(10,B);

    //-----------------------------------------------------------------
    // Генерация оси вращения пульсара и расположения магнитного полюса
    chance_1  = rand()/rand_high_board;
    chance_1 *= 2*pi;

    do {
        chance_2  = rand()/rand_high_board;
        chance_2 *= pi;
        chance_2 -= pi/2.;
        chance_3  = rand()/rand_high_board;
    } while (cos(chance_2) <= chance_3);

    x_axis    = cos(chance_2) * cos(chance_1);
    y_axis    = cos(chance_2) * sin(chance_1);
    z_axis    = sin(chance_2);

    do {
        chance_1  = rand()/rand_high_board;
        chance_1 *= pi;
        chance_1 -= pi/2.;
        chance_3  = rand()/rand_high_board;
    } while (abs(sin(chance_1)) <= chance_3);

    i_incl    = chance_1;

    i_incl = 10./180.*pi;

    //-----------------------------------------------------------------
    // Set the birth time of NS.
    tau = T;
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    // Is it NS? If not then we set mass as BH should have according to
    // work by Fryer, Belczynski, Wiktorowicz at al. (2012) delay collapse

double mass_prog, Z_prog;

    if (proteg->get_c_BAGB()<9.5084) {
        massive = true;
    } else {
        massive = false;
	mass_prog = proteg -> get_mass ();
	Z_prog    = proteg -> get_Z    ();		
	    if (mass_prog > 11 && mass_prog < 30)
		M = 1.1 + 0.2 * exp ((mass_prog - 11.0)/4.) - (2. + Z_prog) * exp(0.4*(mass_prog - 26.0));
	    else if ( mass_prog > 30)
		M = min (33.35 + (4.75 + 1.25*Z_prog)*(mass_prog - 34.), mass_prog - sqrt(Z_prog)*(1.3*mass_prog - 18.35));
	    }

    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Генерируем место расположения сгустков плазмы для старого и
    // молодого пульсаров. А именно 0 - высота в км, 1 - угол в рад.
    // В основном следуем работе Karastergiou, 2008, но с некоторыми
    // видоизменениями

    memset (sparks, 0, sizeof(sparks));

    for (int i = 0; i < 10; i++)	{
        sparks [1][i] = rand() / rand_high_board;
        sparks [1][i] *= 2 * pi;
        sparks [0][i] = rand() / rand_high_board;
        sparks [0][i] *= 100;
        sparks [0][i] += 900;
    }

    int first_group, second_group, third_group;

    first_group  = rand() % 5 + 2.;
    second_group = rand() % 5 + 2.;

    third_group = max(7., 20. - first_group - second_group);


    for (int i = 10; i < first_group+10; i++)	{
        sparks [1][i] = rand() / rand_high_board;
        sparks [1][i] *= 2 * pi;
        sparks [0][i] = rand() / rand_high_board;
        sparks [0][i] *= 300;
    }

    for (int i = 10+first_group; i < first_group+10 + second_group; i++)       {
        sparks [1][i] = rand() / rand_high_board;
        sparks [1][i] *= 2 * pi;
        sparks [0][i] = rand() / rand_high_board;
        sparks [0][i] *= 500;
        sparks [0][i] += 300;
    }

    for (int i = 10 + first_group + second_group; i < 10 + first_group + second_group + third_group; i++)       {
        sparks [1][i] = rand() / rand_high_board;
        sparks [1][i] *= 2 * pi;
        sparks [0][i] = rand() / rand_high_board;
        sparks [0][i] *= 200;
        sparks [0][i] += 800;
    }


   base_mfd = mfd;
   base_lm  = lm; 

}

double NeutronStar::get_M() {
    return M;
}

double NeutronStar::get_R() {
    return R/1e5;
}

//-----------------------------------------------------------------------
// A new constructor basing on data from information file. 
//-----------------------------------------------------------------------
NeutronStar::NeutronStar (double T, MFD * mfd, LM * lm, GD * p_distr, GD * b_distr, double m, double x_ns, double y_ns, double z_ns, double v_x_ns, double v_y_ns, double v_z_ns, double t2) {

tau = T-t2;

x = x_ns;
y = y_ns;
z = z_ns;

v_x = v_x_ns;
v_y = v_y_ns;
v_z = v_z_ns;

M = m;

base_mfd = mfd;
base_lm  = lm; 

    P = p_distr->generate_next();

    B=b_distr->generate_next();
    B=pow(10,B);

    double is_velocity_set;
    double chance_1, chance_2, chance_3;

    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_x += 1e5*lcm/lsec*expon_vel(chance_1);
    } else {
        v_x += -1e5*lcm/lsec*expon_vel(chance_1);
    }

    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_y += 1e5*lcm/lsec*expon_vel(chance_1);
    } else {
        v_y += -1e5*lcm/lsec*expon_vel(chance_1);
    }

    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_z += 1e5*lcm/lsec*expon_vel(chance_1);
    } else {
        v_z += -1e5*lcm/lsec*expon_vel(chance_1);
    }

    //-----------------------------------------------------------------
    // Генерация оси вращения пульсара и расположения магнитного полюса
    chance_1  = rand()/rand_high_board;
    chance_1 *= 2*pi;
     do {
        chance_2  = rand()/rand_high_board;
        chance_2 *= pi;
        chance_2 -= pi/2.;
        chance_3  = rand()/rand_high_board;
    } while (cos(chance_2) <= chance_3);

    //cout<<"phi = "<<chance_1*180./pi<<endl;
    //cout<<"psy = "<<chance_2*180./pi<<endl;

    x_axis    = cos(chance_2) * cos(chance_1);
    y_axis    = cos(chance_2) * sin(chance_1);
    z_axis    = sin(chance_2);

    do {
        chance_1  = rand()/rand_high_board;
        chance_1 *= pi;
        chance_1 -= pi/2.;
        chance_3  = rand()/rand_high_board;
    } while (abs(sin(chance_1)) <= chance_3);

    i_incl    = chance_1;

    i_incl = 10./180.*pi;

    //-----------------------------------------------------------------
    // Генерируем место расположения сгустков плазмы для старого и
    // молодого пульсаров. А именно 0 - высота в км, 1 - угол в рад.
    // В основном следуем работе Karastergiou, 2008, но с некоторыми
    // видоизменениями

    memset (sparks, 0, sizeof(sparks));

    for (int i = 0; i < 10; i++)	{
        sparks [1][i] = rand() / rand_high_board;
        sparks [1][i] *= 2 * pi;
        sparks [0][i] = rand() / rand_high_board;
        sparks [0][i] *= 100;
        sparks [0][i] += 900;
    }

    int first_group, second_group, third_group;

    first_group  = rand() % 5 + 2.;
    second_group = rand() % 5 + 2.;

    third_group = max(7., 20. - first_group - second_group);


    for (int i = 10; i < first_group+10; i++)	{
        sparks [1][i] = rand() / rand_high_board;
        sparks [1][i] *= 2 * pi;
        sparks [0][i] = rand() / rand_high_board;
        sparks [0][i] *= 300;
    }

    for (int i = 10+first_group; i < first_group+10 + second_group; i++)       {
        sparks [1][i] = rand() / rand_high_board;
        sparks [1][i] *= 2 * pi;
        sparks [0][i] = rand() / rand_high_board;
        sparks [0][i] *= 500;
        sparks [0][i] += 300;
    }

    for (int i = 10 + first_group + second_group; i < 10 + first_group + second_group + third_group; i++)       {
        sparks [1][i] = rand() / rand_high_board;
        sparks [1][i] *= 2 * pi;
        sparks [0][i] = rand() / rand_high_board;
        sparks [0][i] *= 200;
        sparks [0][i] += 800;
    }

    if (M<2.8) {
        massive = true;
    } else {
        massive = false;
    }

}
//-----------------------------------------------------------------------
// Check if the pulsar has already crossed the death line
// The check is exactly the same as in Faucher-Giguere & Kaspi (2006) article

bool NeutronStar::is_pulsar_alive(double t) {

    //cout<<get_incl(t)/pi*180<<"\t"<<abs(get_incl(t)/pi*180-90)<<endl;

    //	if (2*pow(get_P(t),11./10.)*pow(get_dot_P(t)/1e-15, -4./10.) >=1.98)
    //		return false;
    //	else if (abs(get_incl(t)/pi*180-90)>0.01)
    //		return true;
    //	else
    //		return false;

    //	if (abs(get_incl(t)/pi*180-90)<=0.1) {
    //		cout<<"It was!"<<endl;
    //		return false;			}



    //return true;

    if (get_B(t)/pow(get_P(t), 2.)>0.12e12) {
        return true;
    } else {
        return false;
    }


}
/*
double neutron_star::get_L(double t) {
double P_res;
P_res  = get_P(t);
return L;
}
*/
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Стандартные для любого класса типа небесного обьекта
// функции, такие как положение объекта, скорость, и метод
// для динамической эволюции положения объекта в Галактике

double NeutronStar::get_position_x() {
    return x;
}

double NeutronStar::get_position_y() {
    return y;
}

double NeutronStar::get_position_z() {
    return z;
}

void NeutronStar::move_to(double T) {
    double result [6];

    T = T - tau;

    result [0] = x;
    result [1] = y;
    result [2] = z;
    result [3] = v_x;
    result [4] = v_y;
    result [5] = v_z;

    //cout<<"A neutron star is moving to "<<T<<endl;

    Runge_Kutta (6, T, &result[0], &diff_equi);

    x   = result [0];
    y   = result [1];
    z   = result [2];
    v_x = result [3];
    v_y = result [4];
    v_z = result [5];

}

double NeutronStar::get_velocity_x (void) {
    return v_x/lcm*lsec/1e5;
}

double NeutronStar::get_velocity_y (void) {
    return v_y/lcm*lsec/1e5;
}

double NeutronStar::get_velocity_z (void) {
    return v_z/lcm*lsec/1e5;
}
//------------------------------------------------------------------

/*
//------------------------------------------------------------------
// Падение магнитного поля в соотвествии со статьёй Pons, 2009
// Формула получена из личной переписки

double neutron_star::get_B (double T) {
double res_B;
T = T - tau;
	if (T<1e6)
		res_B = B * pow(e, -T/tau_ohm) / (1.+alpha*(1.-pow(e, -T/tau_ohm)));
	else
		res_B = B * pow(e, -1e6/tau_ohm) / (1.+alpha*(1.-pow(e, -1e6/tau_ohm)));
return res_B;
}
*/
/*
//------------------------------------------------------------------------
// Метод, позволяющий узнать, виден ли пульсар с Земли. Учитываются такие
// эффекты как биминг, ослабление сигнала с расстоянием, угловое положение
// пульсара на небе. Параметры обзоров пульсаров взяты из Faucher, 2006
//------------------------------------------------------------------------

bool neutron_star::is_pulsar_visible (double t, special_star * sun) {
double chance_1, f, l, b;
double first[3], second[2];
double dist_to_sun;
double lum;

sun->move_to(t);

	if (visible)	{
		dist_to_sun = pow(sun->get_position_x() - x, 2) + pow(sun->get_position_y() - y, 2) + pow(sun->get_position_z() - z, 2);
		lum = get_Lum (t);
		lum /= dist_to_sun;
		// Вектор от солнца к пульсару направленный
		first [0] = sun->get_position_x() - x;
		first [1] = sun->get_position_y() - y;
		first [2] = sun->get_position_z() - z;
		// Вектор от солнца к центру Галактики
		second[0] = - sun->get_position_x();
		second[1] = - sun->get_position_y();

		// Вычисление галактических координат пульсара
		l = (first[0] * second [0] + first[1]*second[1]) / sqrt(pow(first[0], 2) + pow(first[1], 2)) / sqrt(pow(second[0], 2) + pow(second[1], 2));
		l = acos (l)/2./pi*360.;
		b = first[2] / sqrt(pow(first[0], 2) + pow(first[1], 2) + pow(first[2], 2));
		b = asin (b)/2./pi*360.;
*/
//		if (lum>0.5e-3 && abs(b)<15. /*&& l <= 260. && l >= 50*/)
/*		return true;

				}
	else	{
		visible = false;
		return false;	}
return false;
}
*/

/*
//--------------------------------------------------------------------------
// Новая попытка написать функцию описываюшию селекцию.
// Учитываеться сложная форма распределения энергия внутри светового конуса
// а именно нормальное с полушириной зависяшей от напряжённости магнитного
// поля.
//--------------------------------------------------------------------------
double neutron_star::is_pulsar_visible (double t, special_star * sun) {
double res;
double dist_to_sun, w50, lum_0, lum_min = 0.3, omega, theta;
double up_border, down_border, l, b;
double first[3], second[2];
sun->move_to(t);
dist_to_sun = sqrt(pow(sun->get_position_x() - x, 2) + pow(sun->get_position_y() - y, 2) + pow(sun->get_position_z() - z, 2));
                // Вектор от солнца к пульсару направленный
                first [0] = sun->get_position_x() - x;
                first [1] = sun->get_position_y() - y;
                first [2] = sun->get_position_z() - z;


		// Вектор от Солнца к центру Галактики
                second[0] = - sun->get_position_x();
                second[1] = - sun->get_position_y();

                b = first[2] / sqrt(pow(first[0], 2) + pow(first[1], 2) + pow(first[2], 2));
                b = asin (b)/2./pi*360.;


//cout<<get_B(t)<<endl;
		w50 = 1.9*log10(get_B(t))-24.5; //21.5 Было
//cout<<pow(10,w50)<<endl;
		w50 = pow(10, w50) / 180 * pi;

		lum_0 = get_Lum(t) / pow(dist_to_sun, 2);
//		cout<<lum_0<<endl;
		if (lum_0 < lum_min)			{
//	cout<<"Problem with luminusity"<<endl;
			return 0.;
							}
		else 					{
//      cout<<w50<<"\t"<<lum_0<<endl;
		omega = w50 * log(lum_0 / lum_min);
//	cout<<omega<<endl;
		theta = (x_axis * first[0] + y_axis * first[1] + z_axis * first[2])/(dist_to_sun*dist_to_sun);
//	cout<<theta<<"\t"<<cos(omega - i_incl)<<"\t"<<cos(omega + i_incl)<<endl;
		if ((omega + i_incl)>pi/2.)
			down_border = pi/2;
		else
			down_border = omega + i_incl;
		if ((i_incl - omega)<0)
			up_border   = 0;
		else
			up_border   = i_incl - omega;
//	if (t<-160e6)
//cout<<t<<"\t"<<theta<<"\t"<<omega + i_incl<<"\t"<<cos(down_border)<<"\t"<<cos(up_border)<<endl;
		if (((cos(down_border) <= theta && cos(up_border) >= theta) ||
		    (cos(pi - down_border) >= theta && cos(pi - up_border) <= theta))
		    && abs(b)<15.					)	{
//	cout<<"yes"<<endl;
			res = abs(abs(acos (theta)) - abs(i_incl));
			res /= abs(omega);
			res = lum_0 * pow(e, -res);
			return res;								}
		else
			return 0.;			}
return res;
}
*/
bool NeutronStar::is_this_ns (void) {
    return massive;
}

double NeutronStar::get_dist_to_sun(double t, SpecialStar * sun)	{
    double res;
    double x_sun, y_sun, z_sun;

    sun->move_to(t);

    x_sun = sun->get_position_x();
    y_sun = sun->get_position_y();
    z_sun = sun->get_position_z();

    res = sqrt(pow(x - x_sun, 2) + pow(y - y_sun, 2) + pow(z - z_sun, 2));

    return res;
}

double NeutronStar::get_gl(double t, SpecialStar * sun)	{
    double res, l;
    double x_sun, y_sun, z_sun;
    double SP[3], SC[3]; // SP вектор пульср - солнце, SC - вектор солнце - центр Галактики

    sun->move_to(t);

    x_sun = sun->get_position_x();
    y_sun = sun->get_position_y();
    z_sun = sun->get_position_z();

    res = sqrt(pow(x - x_sun, 2) + pow(y - y_sun, 2) + pow(z - z_sun, 2));

    SP [0] = x - x_sun;
    SP [1] = y - y_sun;
    SP [2] = z - z_sun;

    SC [0] = - x_sun;
    SC [1] = - y_sun;
    SC [2] = - z_sun;

    l = acos ((SP[0]*SC[0]+SP[1]*SC[1]) / sqrt(SC[0]*SC[0] + SC[1]*SC[1]) / sqrt(SP[0]*SP[0] + SP[1]*SP[1]));

    if ((SP[0]*sun->get_velocity_x() + SP[1]*sun->get_velocity_y()) < 0) {
        l = 2*pi - l;
    }
    	
    return l*180.0/M_PI;
}

double NeutronStar::get_gb(double t, SpecialStar * sun)	{
    double res, b;
    double x_sun, y_sun, z_sun;
    double SP[3], SC[3]; // SP вектор пульср - солнце, SC - вектор солнце - центр Галактики

    sun->move_to(t);

    x_sun = sun->get_position_x();
    y_sun = sun->get_position_y();
    z_sun = sun->get_position_z();

    res = sqrt(pow(x - x_sun, 2) + pow(y - y_sun, 2) + pow(z - z_sun, 2));

    SP [0] = x - x_sun;
    SP [1] = y - y_sun;
    SP [2] = z - z_sun;

    SC [0] = - x_sun;
    SC [1] = - y_sun;
    SC [2] = - z_sun;
    	
    b = asin (SP[2]/res);

    return b * 180.0 / M_PI;
}

double NeutronStar::get_init_P () {

	return P;

}	

double NeutronStar::get_init_B () {                     

        return B;

}

//---------------------------------------------------------------//
// Получение меры дисперсии для данного наблюения пульсара,
// используется код NE2001 (obsolate)
//---------------------------------------------------------------//

//extern "C" {
//    void dmdsm_ (float *l, float *b, int *ndir, float *dmpsr, float *dist, char *limit, float *sm, float *smtau, float *smtheta, float *smiso);
//}

//---------------------------------------------------------------//
// Получение меры дисперсии для данного наблюения пульсара,
// используется код ymw16 (11-03-2017)
//---------------------------------------------------------------//
extern "C" {
    double d_to_dm (double, double, double);
}

double NeutronStar::get_DM (double t, SpecialStar * sun, float *l, float *b, float *sm) {
    double x_sun, y_sun, z_sun;
    double gl, gb;
    float dist, dmpsr;
    double dmpsr1, dist1;
    //float l, b;
    double SP[3], SC[3]; // SP вектор пульср - солнце, SC - вектор солнце - центр Галактики

    sun->move_to(t);

    x_sun = sun->get_position_x();
    y_sun = sun->get_position_y();
    z_sun = sun->get_position_z();

    dist = sqrt(pow(x - x_sun, 2) + pow(y - y_sun, 2) + pow(z - z_sun, 2));

    SP [0] = x - x_sun;
    SP [1] = y - y_sun;
    SP [2] = z - z_sun;

    SC [0] = - x_sun;
    SC [1] = - y_sun;
    SC [2] = - z_sun;

    *l = acos ((SP[0]*SC[0]+SP[1]*SC[1]) / sqrt(SC[0]*SC[0] + SC[1]*SC[1]) / sqrt(SP[0]*SP[0] + SP[1]*SP[1]));

    if ((SP[0]*sun->get_velocity_x() + SP[1]*sun->get_velocity_y()) < 0) {
        *l = 2*pi - *l;
    }

    *b = asin (SP[2]/dist);

    gl = *l;
    gb = *b;
    dist1 = dist;

    //cout << "TEST: gl, gb, dist: "<< gl << "\t" << gb << "\t"<< dist << endl;
    //cout << "TEST, TEST, TEST: "<< d_to_dm (0.0, 0.0, 500.0) << endl;
    gl = gl*180.0/M_PI;
    gb = gb*180.0/M_PI;
    dist1 = dist * 1000.0;
    
    int ndir = -5;

    float /*sm, */ smtau, smtheta, smiso;

    char limit;
    limit = ' ';

    //cout<<*l<<"\t"<<*b<<"\t"<<dist<<endl;

    //dmdsm_ (l, b, &ndir, &dmpsr, &dist, &limit, sm, &smtau, &smtheta, &smiso);
    //cout << "TEST: gl, gb, dist: "<< gl << "\t" << gb << "\t"<< dist1 << endl;

    dmpsr1 = d_to_dm (gl, gb, dist1);
    //dmpsr1 = 50.0;

    //cout << "Test position: "<<dmpsr1 << "\t" << dmpsr << endl;

    //cout << gl*180.0/M_PI <<"\t"<< gb*180.0/M_PI <<"\t" << dist << "\t" << dmpsr << endl;
    //cout << gl << "\t" << gb << "\t" << dist1 << "\t" << dmpsr1<<endl;


    //exit(0);

    return dmpsr1;
}
