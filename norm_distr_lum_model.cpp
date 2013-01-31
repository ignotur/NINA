#include <cmath>
#include <iostream>
#include "stars.h"

using namespace std;

//-----------------------------------------------------------------------//
// Модель светимости пульсара, при которой энергия в конусе излучения
// распределена по нормальному закону.
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//
// Декларация функций
double S_min (double, double, float, double, double, double, float, TMap *);
//-----------------------------------------------------------------------//

double NeutronStar::is_pulsar_visible (double t, SpecialStar * sun, TMap * T_copy, parametrs_lum * param) {
    double res, time_tmpl;
    float l, b, sm, DM;
    double dist_to_sun, w50, lum_0, lum_min = 0.3, omega, theta;
    double up_border, down_border;
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


    w50 = 1.9*log10(get_B(t))-21.5 + param->get_ds();
    w50 = pow(10, w50) / 180 * pi;

    // Встроенная модель светимости
    time_tmpl = t - tau;
    //	if (time_tmpl > 1000)					{
    lum_0 =  -0.646*log10(time_tmpl) + 7.524 + param->get_dlum();
    lum_0 = pow(10, lum_0);		//		}
    //	else	if (time_tmpl < 1000 && time_tmpl > 0)		{
    //		lum_0 =  -0.646*3. + 7.524 + param->get_dlum();
    //				lum_0 = 1000*param->get_dlum();
    //		lum_0 = pow(10.,lum_0);				}
    //-------------------------------//

    lum_0 /= pow(dist_to_sun, 2);

    DM = get_DM (t, sun, &l, &b, &sm);
    lum_min = S_min (l, b, sm, dist_to_sun, w50, get_P(t), DM, T_copy);


    if (lum_0 < lum_min)			{
        return 0.;
    } else 					{

        // omega - угол, на который мы можем отклониться от центра консуса,
        // чтобы ещё можно было увидеть пульсар.
        omega = w50 * log(lum_0 / lum_min);

        theta = (x_axis * first[0] + y_axis * first[1] + z_axis * first[2])/(dist_to_sun);

        if ((omega + i_incl)>pi/2.) {
            down_border = pi/2;
        } else {
            down_border = omega + i_incl;
        }

        if ((i_incl - omega)<0) {
            up_border   = 0;
        } else {
            up_border   = i_incl - omega;
        }

        if (((cos(down_border) <= theta && cos(up_border) >= theta) ||
                (cos(pi - down_border) >= theta && cos(pi - up_border) <= theta))
                && abs(b)<15.					)	{
            res = abs(abs(acos (theta)) - abs(i_incl));
            res /= abs(omega);
            res = lum_0 * pow(e, -res);
            return res;
        } else {
            return 0.;
        }
    }

    return res;
}


parametrs_lum::parametrs_lum (ifstream * in) {
    *in>>ds;
    *in>>dlum;
}

bool parametrs_lum::is_beam_on(double P) {
    return true;
}

void parametrs_lum::print_description (ostream * out) {
    *out<<"// Используется модель светимости B, а именно модель, где   //"<<endl;
    *out<<"// распределение энергии внутри конуса излучения принято    //"<<endl;
    *out<<"// нормальным                                               //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_lum::print_parametrs   (ostream * out) {
    *out<<"//                  Параметры модели светимости             //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// ds   -  "<<ds<<endl;
    *out<<"// dlum -  "<<dlum<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_lum::print_short      (ostream * out) {
    *out<<"ds   -  "<<ds<<endl;
    *out<<"dlum -  "<<dlum<<endl;
}

double parametrs_lum::get_ds (void) {
    return ds;
}

double parametrs_lum::get_dlum (void) {
    return dlum;
}
