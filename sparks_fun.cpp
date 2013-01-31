#include <cstdlib>
#include <time.h>
#include <cmath>
#include <iostream>
#include "stars.h"
#include <algorithm>

using namespace std;

//----------------------------------------------------------------------//
// Часть класса neutron_star - а именно функция которая работает с
// отдельными сгусткими плазмы, для любового возраста пересчитывает их
// положение в широты и долготы, также функция которая позволяет получить
// профили излучения пульсаров и картину распределения излучающих
// сгустков плазмы в магнетосфере в проекции ортогональной конусу
// излучения, а также функция возращающия светимость пульсара
//
//  Автор: Igoshev Andrey
//  Научный руководитель: Холтыгин А.Ф.
//  e-mail: igoshev-andrei@rambler.ru
//  Написание начато: 10.03.2011
//-----------------------------------------------------------------------//

//-----------------------------------------------------------------------//
// Декларация функций
double S_min (double, double, float, double, double, double, float, T_map *);
//-----------------------------------------------------------------------//

//-----------------------------------------------------------------------//
// Функция выводящия расположение сгустков плазмы - преимущественно для
// отладки программы

void neutron_star::show_pos_sparks (void) {

    for (int i = 10; i < 30; i++)	{
        cout<<sparks[1][i]<<"\t"<</*i<<") rho: "<<*/sqrt(9*pi*sparks[0][i]*1.e5/(2.*light_velocity*get_P(tau)))<<"\t"/*<<", s: "<<sparks[1][i]*/<<endl;
    }

}
//-----------------------------------------------------------------------//

//-----------------------------------------------------------------------//
// Функция считающия угол между ближайшим по широте сгустком плазмы и
// линией зрения, также переводит координаты сгустков в широты
// Возвращает отношение этого минимального угла к полуширине распределения
// Дополнительно возвращает полуширину пика, через третий аргумент
// по статье Karastergiou, 2008

double neutron_star::get_angle_btw_l_of_s_sparks (double t, special_star * sun, double * w, parametrs_lum * param) {
    int n_begin;
    int n_end;
    int j_tmpl;
    double latitudes [20];
    double rho       [20];

    if (get_P (t) < 0.15)	{
        n_begin = 0;
        n_end  = 9;
    } else			{
        n_begin = 10;
        n_end   = 29;
    }

    for (int i = n_begin; i < n_end; i++)	{
        rho [i - n_begin] = sqrt(9*pi*sparks[0][i]*1.e5/(2.*light_velocity*get_P(t)))*pi/180./2.;
        latitudes [i - n_begin] = cos(rho[i - n_begin]) / (sqrt(1. - pow(sin(sparks[1][i]) * sin(rho[i-n_begin]), 2)));
        latitudes [i - n_begin] = acos (latitudes [i - n_begin]);
        latitudes [i - n_begin] += i_incl;
    }

    sun->move_to(t);

    double dist_to_sun;
    double first[3], second[2];
    double theta;
    double compr;

    dist_to_sun = sqrt(pow(sun->get_position_x() - x, 2) + pow(sun->get_position_y() - y, 2) + pow(sun->get_position_z() - z, 2));

    // Вектор от солнца к пульсару направленный
    first [0] = sun->get_position_x() - x;
    first [1] = sun->get_position_y() - y;
    first [2] = sun->get_position_z() - z;

    theta = (x_axis * first[0] + y_axis * first[1] + z_axis * first[2])/(dist_to_sun*dist_to_sun);
    theta = acos(theta);
    theta = pi/2. - theta;

    compr = pi;

    for (int i = n_begin; i < n_end; i++) 	{
        if  (min(compr, abs(latitudes[i - n_begin] - theta)) < compr ||  min(compr, abs(-latitudes[i - n_begin] - theta)) > compr) {
            compr = min(compr, abs(latitudes[i - n_begin] - theta));
        }

        compr = min(compr, abs(-latitudes[i - n_begin] - theta));
        j_tmpl= i;
    }

    compr =  - compr/0.042761/param->get_ds()/sqrt(sparks[0][j_tmpl]/10./get_P(t));
    *w     = 0.042761*param->get_ds()*sqrt(sparks[0][j_tmpl]/10./get_P(t));

    return compr;
}

//----------------------------------------------------------------------//

//----------------------------------------------------------------------//
// Функция выводящия значения отчётов для среднего профиля пульсара
// Код по большей части является повторение функции get_angle_btw_l_of_s_sparks
//----------------------------------------------------------------------//
void neutron_star::show_pulse_profile (double t, special_star * sun, parametrs_lum * param) {
    int n_begin;
    int n_end;
    double latitudes [20];
    double longtitude[20];
    double rho       [20];

    if (get_P (t) < 0.15)   {
        n_begin = 0;
        n_end  = 9;
    } else                    {
        n_begin = 10;
        n_end   = 29;
    }

    for (int i = n_begin; i < n_end; i++)   {
        rho [i - n_begin] = sqrt(9*pi*sparks[0][i]*1.e5/(2.*light_velocity*get_P(t)))*pi/180./2.;
        latitudes [i - n_begin] = cos(rho[i - n_begin]) / (sqrt(1. - pow(sin(sparks[1][i]) * sin(rho[i-n_begin]), 2)));
        latitudes [i - n_begin] = acos (latitudes [i - n_begin]);
        latitudes [i - n_begin] += i_incl;
        longtitude[i - n_begin] = sin(sparks[1][i]) * sin(rho[i]);
        longtitude[i - n_begin] = asin(longtitude[i - n_begin]);
    }

    sun->move_to(t);

    double dist_to_sun;
    double first[3], second[2];
    double theta;
    double compr;

    dist_to_sun = sqrt(pow(sun->get_position_x() - x, 2) + pow(sun->get_position_y() - y, 2) + pow(sun->get_position_z() - z, 2));

    // Вектор от солнца к пульсару направленный
    first [0] = sun->get_position_x() - x;
    first [1] = sun->get_position_y() - y;
    first [2] = sun->get_position_z() - z;

    theta = (x_axis * first[0] + y_axis * first[1] + z_axis * first[2])/(dist_to_sun*dist_to_sun);
    theta = acos(theta);
    theta = pi/2. - theta;

    double one_point;
    double azim_cent = pi/2.;

    for (int i = 0; i < 1800; i++)	{
        one_point = 0;

        for (int j = n_begin; j < n_end; j++) {
            one_point += pow(e, -acos(cos(pi/2.+longtitude[j-n_begin]-i*1.74533e-3) * cos(abs(theta - latitudes[j - n_begin])))/0.042761/param->get_ds()/sqrt(sparks[0][j]/10./get_P(t)));
        }

        cout<<i*0.1<<"\t"<<one_point<<endl;
    }
}


//------------------------------------------------------------------------//
// Функция светимости, говорит виден ли пульсар. Если не виден возвращает
// ноль, если виден возвращает его яркость на длине волны 1400 Мгц в яньских
// так же учитывается галактические координаты пульсара
// По статье Fan, 2001
//------------------------------------------------------------------------//
double neutron_star::is_pulsar_visible (double t, special_star * sun, T_map * T_copy, parametrs_lum * param) {
    double res, w;
    float  l, b, sm, DM;
    double dist_to_sun, lum_0, lum_min = 0.3; // Пока оценочное значение
    double rho_6 = 0.5; // Выбрано из описаного в статье интервала
    double first[3], second[2];
    double t_tmp;
    t_tmp = t;

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
    b = asin (b)/pi*180.;

    //		res = 6.25e14*delta_lum*2.*pi*pi/light_velocity*pow(rho_6, 8./7.)*pow(get_P(t)*get_B(t)/1.e12, -2./7.);

    res = 4e5*param->get_dlum();

    res /= pow(dist_to_sun, 2);
    res *= pow(e, get_angle_btw_l_of_s_sparks (t_tmp, sun, &w, param));

    // Вызов функций, которые дадут в итоге минимальную чувствительность
    // радиотелескопа
    sun->move_to(-t);

    w = w/2./pi*get_P(t);
    DM = get_DM (t, sun, &l, &b, &sm);

    //cout<<"W      - "<<"\t"<<w<<"\t"<<get_P(t)<<endl;

    lum_min = S_min (l, b, sm, dist_to_sun, w, get_P(t), DM, T_copy);

    //cout<<"lum_min - "<<"\t"<<lum_min<<endl;

    if (res >=lum_min && abs(b)<=15.0) {
        return res;
    } else {
        return 0;
    }

    //neutron_star::get_angle_btw_l_of_s_sparks(double&, special_star**)
    //neutron_star::get_angle_btw_l_of_s_sparks(double, special_star*)


}


parametrs_lum::parametrs_lum (ifstream * in) {
    *in>>ds;
    *in>>dlum;
}

void parametrs_lum::print_description (ostream * out) {
    *out<<"// Используется модель светимости А, а именно модель со сгу-//"<<endl;
    *out<<"// стками плазмы, расположенными на различных высотах над   //"<<endl;
    *out<<"// нейтронной звездой                                       //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_lum::print_parametrs   (ostream * out) {
    *out<<"//                  Параметры модели светимости             //"<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
    *out<<"// ds   -  "<<ds<<endl;
    *out<<"// dlum -  "<<dlum<<endl;
    *out<<"//----------------------------------------------------------//"<<endl;
}

void parametrs_lum::print_short   (ostream * out) {
    *out<<"ds   -  "<<ds<<endl;
    *out<<"dlum -  "<<dlum<<endl;
}

double parametrs_lum::get_ds (void) {
    return ds;
}

double parametrs_lum::get_dlum (void) {
    return dlum;
}
