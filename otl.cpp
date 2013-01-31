#include <iostream>
#include "stars.h"
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstring>

double phi (double, double, double);

using namespace std;
int main () {

    ifstream input     ("input.txt");

    parametrs_B   param_B   (&input);
    parametrs_lum param_lum (&input);

    double T = -100e6;              // начало
    int number_stars = 7;		// темп звездообразования (звёзд в тысячалетие)
    int number_millenium = -T/1e3;  // количество тысячалетий

    cout<<"//----------------------------------------------------------//"<<endl;
    cout<<"// Пульсарная популяция. Версия 0.80см соглашение о версиях //"<<endl;
    cout<<"// Автор: Игошев Андрей, научный руководитель: А.Ф. Холтыгин//"<<endl;
    cout<<"// e-mail: igoshev-andrei@rambler.ru СПбГУ, 2010-2011       //"<<endl;
    cout<<"//----------------------------------------------------------//"<<endl;
    param_B.print_description  (&cout);
    param_B.print_parametrs    (&cout);
    param_lum.print_description(&cout);
    param_lum.print_parametrs  (&cout);

    cout<<"// Параметры популяциии:                                          "<<endl;
    cout<<"// T_start "<<T<<endl;
    cout<<"// star formation rate "<<number_stars<<endl;
    cout<<"//----------------------------------------------------------//"<<endl;
    cout<<"// Инициализация расчётов.                                  //"<<endl;
    srand(time(0));
    T_map T_copy;

    cout<<"// Инициализация расчётов закончена.                        //"<<endl;
    cout<<"//----------------------------------------------------------//"<<endl;

    ofstream otl_rt ("otl_rt.txt");
    ofstream otl_P  ("otl_P.txt");
    ofstream otl_dot_P ("otl_dot_P.txt");
    ofstream otl_B  ("otl_B.txt");
    ofstream otl_pos ("otl_pos.txt");

    double tmpl, initial;
    special_star sun;
    double now = 0, shift;
    double P, dot_P, B, x, y, z;
    double dist_to_sun, lumin;
    int counter = 0;
    int rand_shift;                 // случайное время рождения внутри тысячелетия
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime(&rawtime);


    star_OB      * ancester;
    neutron_star * descendant;

    double x_i, y_i, z_i;
    double v_x, v_y, v_z, incl;


    T = -100e6;
    //        sun.move_to(T);
    //	for (int i = 0; i < 5000; i++)	{
    //	ancester = new star_OB (T, &sun);
    //		x_i = ancester->get_position_x();
    //		y_i = ancester->get_position_y();
    //		z_i = ancester->get_position_z();
    //
    //	otl_pos<<x_i<<"\t"<<y_i<<"\t"<<z_i<<endl;

    //	delete ancester;	}

    ancester = new star_OB (T, &sun);

    //cout<<"Look here "<<param_lum.is_beam_on(10.)<<endl;

    descendant = new neutron_star(T, ancester, &param_B);
    //	ancester->move_to(T);

    for (int i = 30; i < 80; i++)	{

        //	rand_shift = rand()%500;
        //	ancester = new star_OB (T + rand_shift, &sun);

        //	descendant = new neutron_star(T, ancester, &param_B);

        //cout<<"Look here "<<param_lum.is_beam_on(10.)<<endl;

        //	otl<<log10(descendant->get_B(T))<<endl;
        //	otl<<descendant->get_R()<<endl;

        //	x_i = descendant->get_velocity_x();
        //	y_i = descendant->get_velocity_y();
        //	z_i = descendant->get_velocity_z();

        //	x_i = sqrt(pow(x_i,2)+pow(y_i,2)+pow(z_i,2));

        //	otl_pos<<x_i<<endl;
        //	cout<<"+1"<<endl;
        //	delete ancester;
        //	delete descendant;

        //	ancester->move_to (pow(10, 2*0.09*i-1.));
        //	ancester->move_to (-pow(10, 2*0.09*i-1.));

        //	x = sun.get_position_x();
        //	y = sun.get_position_y();
        //	z = sun.get_position_z();

        //	v_x = sun.get_velocity_x();
        //	v_y = sun.get_velocity_y();
        //	v_z = sun.get_velocity_z();

        //	x = sun.get_position_x();
        //	y = sun.get_position_y();

        //        sun.move_to(T+1.8e6*i);
        //	z = ancester->get_position_z();

        //	otl<<pow(10, 2*0.09*i-1.)<<"\t"<<sqrt(pow(x-x_i, 2) + pow(y-y_i, 2) + pow (z-z_i, 2))<<endl;

        //	initial =  pow(ancester->get_velocity_x(), 2) + pow(ancester->get_velocity_y(), 2) + pow(ancester->get_velocity_z(), 2) + 2*phi(ancester->get_position_x(), ancester->get_position_y(), ancester->get_position_z());

        //	otl<<log10(ancester->get_Flux())<<endl;
        //	otl<<ancester->get_mass()<<"\t"<<ancester->get_R(now)<<endl;

        //	for (int i = 0; i < 1000; i++)	{

        //	descendant->move_to (T+1e3*i);

        //	otl_pos<<1e3*i<<"\t"<<descendant->get_position_x()<<"\t"<<descendant->get_position_y()<<"\t"<<descendant->get_position_z()<<endl;
        //	otl_pos<<sqrt(x*x+y*y)<<"\t"<<sqrt(v_x*v_x+v_y*v_y)/1e5/lcm*lsec<<endl;


        //---------------------------------------------------------------------------------------------------------//
        // B, P, incl, tau
        B    = descendant->get_B   (T+pow(10, 0.09*i-1));
        P    = descendant->get_P   (T+pow(10, 0.09*i-1));
        incl = descendant->get_incl(T+pow(10, 0.09*i-1));
        //		if (descendant->is_pulsar_alive(T+pow(10, 0.09*i-1)))
        //		if (descendant->is_pulsar_alive(T+pow(10, 0.09*i-1)))
        otl_P<<pow(10, 0.09*i-1)<<"\t"<<B<<"\t"<<P<<"\t"<<"\t"<<P/2./descendant->get_dot_P(T+pow(10, 0.09*i-1))/3.2e7/*<<"\t"<<P*P/2./(1.3e-39*B*B*cos(incl)*cos(incl))/3.2e7*/<<"\t"<<descendant->get_dot_P(T+pow(10, 0.09*i-1))/*<<"\t"<<B*B*cos(incl)*cos(incl)/P/3.2e19/3.2e19*/<<endl;
        //---------------------------------------------------------------------------------------------------------//



        //	dot_P = descendant->get_dot_P(T+pow(10, 0.09*i-1));

        //        otl_pos<<1.8e6*i<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<v_x<<"\t"<<v_y<<"\t"<<v_z<<endl;

        //	otl_rt<<P/2./dot_P/3600./24./365.24<<"\t"<<pow(10, 0.09*i-1)<<endl;
        //	otl_B<<P/2./dot_P/3600./24./365.24<<"\t"<<B<<endl;
        //	otl_P<<P/2./dot_P/3600./24./365.24<<"\t"<<P<<endl;
        //	otl_dot_P<<P/2./dot_P/3600./24./365.24<<"\t"<<dot_P<<endl;
        //	ancester->move_to (pow(10, 2*0.09*i-1.));
        //	tmpl = (pow(ancester->get_velocity_x(), 2) + pow(ancester->get_velocity_y(), 2) + pow(ancester->get_velocity_z(), 2) + 2*phi(ancester->get_position_x(), ancester->get_position_y(), ancester->get_position_z())-initial)/initial;
        //	otl<<pow(10, 2*0.09*i-1.)<<"\t"<<tmpl<<endl;

        //	delete descendant;
        //	delete ancester;
    }

    delete descendant;
    delete ancester;



    //	shift  = ancester->get_time_on_MS();
    //	shift += ancester->get_t_He      ();

    //	ancester->move_to (T + shift + rand_shift);
    //	descendant = new neutron_star (T + shift + rand_shift, ancester, &param_B);
    //	cout<<"Has it alived? "<<descendant->is_pulsar_alive(now)<<endl;

    //	descendant->move_to(now);
    //	lumin = descendant->is_pulsar_visible(now, &sun, &T_copy, &param_lum);
    //	cout<<"lumin - "<<lumin<<endl;
    //	cout<<"Binit - "<<descendant->get_B(T + shift + rand_shift)<<endl;
    //	cout<<"Binit - "<<descendant->get_B(T + shift + rand_shift + param_B.get_time_II())<<endl;
    //	for (int i = 0; i < 1000; i++)
    //	otl<<pow(10, 0.009*i-1.)<<"\t"<<descendant->get_B(pow(10, 0.009*i-1.)+ T + shift + rand_shift)<<endl;
    //	otl<<pow(10, 0.009*i-1.)<<"\t"<<descendant->get_P(pow(10, 0.009*i-1.)+ T + shift + rand_shift)<<"\t"<<descendant->get_dot_P(pow(10, 0.009*i-1.) + T + shift + rand_shift)<<endl;
    //	delete descendant;
    //delete ancester;

    return 0;
}
