#ifndef STAR
#define STAR

#include "lum_model.h"
#include "field_decay_model.h"

double const light_velocity = 2.9979250e10;  // Скорость света в см/сек
double const pi      = 3.1415926;
double const rand_high_board = 2.14665e+9;
double const e       = 2.718281828;
double const M_sol   = 1.989e33;                // Масса солнца в граммах
double const R_sol   = 6.9599e10;               // Радиус солнца в сантиметрах
double const G_cgs   = 6.67e-8;                 // Гравитационная постоянная в см^3 /(г*сек^2)
double const G       = 2.2608e-57;              // Гравитационная постоянная в кпк^3/(г*год^2)
double const lcm     = 3.2407789e-22;           // 1 см  в кпк,   для внутренних переводов
double const lsec    = 3.1688955e-8;            // 1 сек в годах, для внутренних переводов

//------------------------------------------------------------//
// Константы по которым будет производиться градиентный спуск
//extern double tau_ohm/*   = 2.3e5*/;                   // характерное время убывания поля на поверхности нейтронной звезды
//extern double alpha/*     = 5*/;                       // коэффициент затухания магнитного поля
//extern double delta_s/*   = 0.1*/;                     // Коэффициент определяющий угловые размеры сгустка плазмы
//extern double delta_lum/* = 1.*/;                      // Коэффициент определяющий светимость одного сгустка плазмы
//------------------------------------------------------------//

//--------------------------------------------------//
// Класс служит для хранения карты распределения
// температур по небу на частоте 1.4 Ггц
// Используется для ускорения расчётов.
//--------------------------------------------------//

class T_map {
public:
    int size;
    double Tb [1038961];
    T_map ();
    double get_Tb (int);
private:

};

//--------------------------------------------//
// Class which contains P and B initial parame-
// ters  for distribution.
//--------------------------------------------//


class P_distr {
public:
    P_distr (ifstream *);
    void print_param (ostream * );
    double a();
    double b();
private:
    double ain;
    double bin;
};

class B_distr {
public:
    B_distr (ifstream *);
    void print_param (ostream *);
    double a();
    double b();
private:
    double ain;
    double bin;
};

//----------------------------------------------//
// Класс для расчёта положения солнца
//----------------------------------------------//

class special_star {
public:
    double x, y, z, v_x, v_y, v_z, t; // эпоха для которой посчитано положение
    int use;
    ofstream out_err;
    special_star();
    //~special_star();
    double get_position_x();
    double get_position_y();
    double get_position_z();
    double get_velocity_x();
    double get_velocity_y();
    double get_velocity_z();
    double get_theta     ();
    void move_to(double);
private:
};
// Описание переменных
// M - масса в массах солнца в граммах
// x, y, z - координаты звезды в кпк
// v_x, v_y, v_z - скорость звезды в кпк/год
// F - магнитный поток на поверхности звезды в Гс см^2
// Z - металличность
// ksi - log10(Z/0.02) - отношение металличности к солнечной
// L - момент импульса находяшийся в звезде в г см^2 / сек^2
class star_OB {
public:
    double M, R, x, y, z, v_x, v_y, v_z, F, P, Z, ksi, L, t_bgb, t_ms, t_inf, t_HeI, t_He, M_c_SN;
    double M_c_DU, M_c_BGB, t_BGB, M_c_HeF, M_c_HeI, M_c_BAGB;
    double tau; // Время рождения
    star_OB  (double, special_star*);
    //~star_OB ();
    double get_position_x();
    double get_position_y();
    double get_position_z();
    double get_velocity_x();
    double get_velocity_y();
    double get_velocity_z();
    void   move_to (double);
    double get_Flux      ();
    double get_mass      ();
    double get_time_on_MS();
    double get_R   (double);
    double get_Z         ();
    double get_L         ();
    double get_t_bgb     ();
    double get_t_inf     ();
    double get_t_HeI     ();
    double get_M_c_SN    ();
    double get_t_He      ();
    double get_c_HeI     ();
    double get_c_DU      ();
    double get_c_BGB     ();
    double get_c_HeF     ();
    double get_c_BAGB    ();
protected:
    void   set_position(double, double, double);
    void   set_velocity(double, double, double);
private:

};
/*
class close_double_star_OB {
public:
star_OB *first, *second;
double x, y, z, v_z, v_y, v_z;
double a, e;
private:
}
*/
class neutron_star {
public:
    parametrs_B * paramet_B;
    double M, R, x, y, z, v_x, v_y, v_z, B, P, L;
    double tau; // Время рождения
    bool visible, massive;
    double sparks [2][30];
    double i_incl, x_axis, y_axis, z_axis;
    neutron_star  (double, star_OB *, parametrs_B *, P_distr*, B_distr*);
    //~neutron_star ();
    double get_M          ();
    double get_R          ();
    double get_P   	(double/*, parametrs_B **/);
    double get_dot_P(double/*, parametrs_B **/);
    //double get_L    (double);
    double get_position_x ();
    double get_position_y ();
    double get_position_z ();
    double get_velocity_x ();
    double get_velocity_y ();
    double get_velocity_z ();
    double get_dist_to_sun(double, special_star *);
    void   move_to  (double);
    double get_B    (double/*, parametrs_B **/);
    double get_incl (double);
    //double get_Lum  (double);
    //double get_nu_0 (double);	// Циклическая частота вращения пульсара sec^-1
    //double get_nu_1 (double);	// Первая производная циклической частоты вращения пульсара sec^-2
    //double get_nu_2 (double);	// Вторая производная циклической частоты вращения пульсара sec^-3
    void show_pos_sparks (void);
    void show_pulse_profile (double t, special_star *, parametrs_lum *);
    bool is_pulsar_alive   (double);
    bool is_this_ns();		// Проверка, какая масса у родившейся нейтронной звезды
    double is_pulsar_visible           (double, special_star * , T_map *, parametrs_lum *);
    double get_angle_btw_l_of_s_sparks (double, special_star *, double *, parametrs_lum *);
    double get_DM (double, special_star *, float *, float *, float *);
    //double T_sky  (double, double);
    //double S_min  (double, double, float, double, double, double, float);
private:

};

#endif
