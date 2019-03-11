#include <cmath>
#include <iostream>
#include <cstdlib>
#include "stars.h"

using namespace std;

//-----------------------------------------------------------------------//
// Модель светимости пульсара, при которой энергия в конусе излучения
// распределена равномерно, как это было сделано у Faucher ...
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
// Декларация функций
double S_min (double, double, float, double, double, double, float, TMap *);
double norm_distr();
//-----------------------------------------------------------------------//

double LMFlat::is_pulsar_visible (double t, SpecialStar * sun, TMap * T_copy, double x, double y, double z, double i_incl, double P, double dot_P, float DM) {
    float l, b, sm;
    double L, eps_P = -1.5, eps_dot_P = 0.5, L_corr=0.8, L_0 = 0.5e-2;
    double dist_to_sun, lum_min, w50;
    double first[3], second[2];

//    sun->move_to(t);
    dist_to_sun = sqrt(pow(sun->get_position_x() - x, 2) + pow(sun->get_position_y() - y, 2) + pow(sun->get_position_z() - z, 2));
    //sun->move_to(-t);

    //                w50 = 0.05*get_P(t);

    w50 = 6.81e-3*sqrt(P)/sin(i_incl);

    //		DM = get_DM (t, sun, &l, &b, &sm);
    //cout<<"Same pulsar"<<endl;
    //cout<<l<<"\t"<<b<<endl;


    first [0] = x - sun->get_position_x();
    first [1] = y - sun->get_position_y();
    first [2] = z - sun->get_position_z();


    // Вектор от Солнца к центру Галактики
    second[0] = - sun->get_position_x();
    second[1] = - sun->get_position_y();

    b = first[2] / dist_to_sun;
    b = asin (b)/2./pi*360.;

    l = atan2(second[0]*first[1]-second[1]*first[0], first[0]*second[0]+first[1]*second[1])*180./3.1415926;

    if (l<0) {
        l=360.+l;
    }

    L=log(L_0*pow(P, eps_P)*pow(dot_P/1e-15, eps_dot_P))+L_corr*norm_distr();
    L = pow(2.71828, L)/pow(dist_to_sun, 2.);


    if (dist_to_sun < 25 && abs(b) < 15 && (l<=50 || l>=230) && L>=5.e-6)  	{
        //	DM=15*dist_to_sun;
//        DM = get_DM (t, sun, &l, &b, &sm);
        lum_min = S_min (l, b, sm, dist_to_sun, w50, P, DM, T_copy);
	//cout << "Actual luminosity is -- "<<L<<"\t , "<<dist_to_sun <<endl;
    } else {
        lum_min = 1e9;
    }




    //cout<<l/180.*3.1415926<<"\t"<<b/180.*3.1415926<<endl;
    //cout<<"Finita la comedia"<<endl;

    //cout<<L<<"\t"<<lum_min<<endl;
    if (L > lum_min /*&& abs(b)<15. && (l>=230. || l<=50.)*/) {
        return L;
    } else {
        return 0;
    }
}

//--------------------------------------------//
// Do not use this function twice!!!
// For singular use only!
//--------------------------------------------//
bool LMFlat::is_beam_on(double P) {
    double f, chance_1;
    f=0.09*pow(log10(P/10.0),2.)+0.03;
    chance_1 = rand () / rand_high_board;

    if (f >= chance_1) {
        return 1;
    } else {
        return 0;
    }

}


LMFlat::LMFlat (vector <double> *) {
}

void LMFlat::print_description (ostream * out) {
    *out<<"#// Model C is the model where energy distribution in the    //"<<endl;
    *out<<"#// cone is flat i.e. same as in Manchester et al. (2006)    //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

void LMFlat::print_parameters   (ostream * out) {
    *out<<"#//            This model has no parameters.                 //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

