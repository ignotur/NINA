#ifndef STAR
#define STAR


#include <vector>
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

void print_head     (ostream *);
void print_help     ();
void print_exclusion();
void print_param (ostream *);
void print_error_no_RDF(char);
void print_error_parameters_not_enough ();
void print_error_flag_non_recognised (char val);
void print_error_no_MFD (char);
void print_error_no_LM  (char);

double version ();

/*struct list {
list * next;
list * previous;
double value1;
double value2;
int num;
};*/

void input_syntax (ifstream *,  char *, bool *, char *, std::vector<double>*, 
                   char *, std::vector<double>*, char *, std::vector<double>*, char *,
		   std::vector<double>*, double *, double *);

void print_param (ostream *,  char, bool, char, vector<double>, 
		   char, vector<double>, char, vector <double>, char,
		   vector<double>, double, double);

//--------------------------------------------------//
// Класс служит для хранения карты распределения
// температур по небу на частоте 1.4 Ггц
// Используется для ускорения расчётов.
//--------------------------------------------------//

class TMap {
public:
    int size;
    double Tb [1038961];
    TMap ();
    double get_Tb (int);
private:

};

//----------------------------------------------//
// Класс для расчёта положения солнца
//----------------------------------------------//

class SpecialStar {
public:
    double x, y, z, v_x, v_y, v_z, t; // эпоха для которой посчитано положение
    int use;
    ofstream out_err;
    SpecialStar();
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

class RDF;

// Описание переменных
// M - масса в массах солнца в граммах
// x, y, z - координаты звезды в кпк
// v_x, v_y, v_z - скорость звезды в кпк/год
// F - магнитный поток на поверхности звезды в Гс см^2
// Z - металличность
// ksi - log10(Z/0.02) - отношение металличности к солнечной
// L - момент импульса находяшийся в звезде в г см^2 / сек^2
class OBStar {
public:
    double M, R, x, y, z, v_x, v_y, v_z, F, P, Z, ksi, L, t_bgb, t_ms, t_inf, t_HeI, t_He, M_c_SN;
    double M_c_DU, M_c_BGB, t_BGB, M_c_HeF, M_c_HeI, M_c_BAGB;
    double tau; // Время рождения
    OBStar  (double, SpecialStar*, RDF *, bool);
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
//--------------------------------------------------//
// Parent class for generation of distribution 
//--------------------------------------------------//
class GD {
public:
virtual double generate_next () {};
virtual void print_param (ostream *){};
};
//--------------------------------------------------//
// Child class for generation of gaussian distribution 
//--------------------------------------------------//
class GDGauss : public GD {
private:
vector <double> * values;
public:
GDGauss (vector <double>*);
double generate_next ();
void   print_param (ostream *);
};

//--------------------------------------------------//
// Child class for generation of multiple gaussian distribution 
//--------------------------------------------------//

class GDMGauss : public GD {
private:
vector <double> * values;
public:
GDMGauss (vector <double>*);
double generate_next ();
void print_param (ostream *);
};
//--------------------------------------------------//
// Parent class for model of magnetic field decay 
//--------------------------------------------------//
class MFD {
public:
virtual double get_P     (double, double, double, double) {};
virtual double get_dot_P (double, double, double, double) {};
virtual double get_B     (double, double, double, double) {};
virtual double get_incl  (double, double, double, double) {};
virtual void print_description (ostream *) {};
virtual void print_parameters  (ostream *) {};
};
//--------------------------------------------------//
// Child class for constant magnetic field 
//--------------------------------------------------//

class MFDConst : public MFD {
public:
MFDConst (vector <double> *);
double get_P     (double, double, double, double);
double get_dot_P (double, double, double, double);
double get_B     (double, double, double, double); 
double get_incl  (double, double, double, double);
void print_description (ostream *);
void print_parameters  (ostream *);
};
//--------------------------------------------------//
// Child class for step magnetic field decay
//--------------------------------------------------//

class MFDStep : public MFD {
public:
MFDStep (vector <double> *);
double get_P     (double, double, double, double);
double get_dot_P (double, double, double, double);
double get_B     (double, double, double, double); 
double get_incl  (double, double, double, double);
void print_description (ostream *);
void print_parameters  (ostream *);
private:
double t_1, t_2, step;
};
//--------------------------------------------------//
// Child class for Pons-like magnetic field decay
//--------------------------------------------------//

class MFDPons : public MFD {
public:
MFDPons (vector <double> *);
double get_P     (double, double, double, double);
double get_dot_P (double, double, double, double);
double get_B     (double, double, double, double); 
double get_incl  (double, double, double, double);
void print_description (ostream *);
void print_parameters  (ostream *);
private:
double tau_ohm, alpha;
};

//--------------------------------------------------//
// Child class for exponential magnetic field decay
//--------------------------------------------------//

class MFDExpon : public MFD {
public:
MFDExpon (vector <double> *);
double get_P     (double, double, double, double);
double get_dot_P (double, double, double, double);
double get_B     (double, double, double, double); 
double get_incl  (double, double, double, double);
void print_description (ostream *);
void print_parameters  (ostream *);
private:
double tau_ohm;
};

//--------------------------------------------------//
// Child class for old Pons' magnetic field decay
//--------------------------------------------------//

class MFDOldPons : public MFD {
public:
MFDOldPons (vector <double> *);
double get_P (double, double, double, double);
double get_dot_P (double, double, double, double);
double get_B (double, double, double, double); 
double get_incl (double, double, double, double); 
void print_description (ostream *);
void print_parameters (ostream *); 
private:
double tau_ohm, tau_hall;
double * pointer_b, * pointer_delta;
};


//--------------------------------------------------//
// Parent class for model of luminosity 
//--------------------------------------------------//
class LM {
public:
virtual double is_pulsar_visible (double, SpecialStar *, TMap *, double, double, double, double, double, double, float) {};
virtual bool   is_beam_on (double) {};
virtual void print_description (ostream *) {};
virtual void print_parameters  (ostream *) {}; 
};

//--------------------------------------------------//
// Child class for model of luminosity (flat distribution)
//--------------------------------------------------//
class LMFlat : public LM {
public:
double is_pulsar_visible (double, SpecialStar *, TMap *, double, double, double, double, double, double, float);
bool   is_beam_on (double);
LMFlat (vector <double> *);
void print_description (ostream *);
void print_parameters  (ostream *); 
private:
double eps_P, eps_dot_P, L_0;
};

//--------------------------------------------------//
// Child class for model of luminosity (exponential distribution)
//--------------------------------------------------//
class LMExpon : public LM {
public:
double is_pulsar_visible (double, SpecialStar *, TMap *, double, double, double, double, double, double, float);
bool   is_beam_on (double);
LMExpon (vector <double> *);
void print_description (ostream *);
void print_parameters  (ostream *); 
private:
double ds, dlum, x_axis, y_axis, z_axis;
bool   does_axis_set;
};

//--------------------------------------------------//
// Parent class for radial distribution of stars 
//--------------------------------------------------//
class RDF {
public:
virtual double rho (double) {};
virtual void   print_description (ostream *){};
};

//--------------------------------------------------//
// Child class for radial distribution of stars (Faucher) model A
//--------------------------------------------------//
class RDFFaucher : public RDF {
public:
double rho (double);
void   print_description (ostream *);
};

//--------------------------------------------------//
// Child class for radial distribution of stars (Kruit) model B
//--------------------------------------------------//
class RDFKruit : public RDF {
public:
double rho (double);
void   print_description (ostream *);
};

//--------------------------------------------------//
// Child class for radial distribution of stars (B0) model C
//--------------------------------------------------//
class RDFB0 : public RDF {
public:
double rho (double);
void   print_description (ostream *);
};

//--------------------------------------------------//
// Child class for radial distribution of stars (SN) model D
//--------------------------------------------------//
class RDFSN : public RDF {
public:
double rho (double);
void   print_description (ostream *);
};

//--------------------------------------------------//
// Child class for radial distribution of stars (Pulsars old) model E
//--------------------------------------------------//
class RDFPuls : public RDF {
public:
double rho (double);
void   print_description (ostream *);
};

//--------------------------------------------//
// Class which contains P and B initial parame-
// ters  for distribution.
//--------------------------------------------//


class PDistr {
public:
    PDistr ();
    void Set (ifstream *);
    void Set (double, double);
    void print_param (ostream * );
    double a();
    double b();
private:
    double ain;
    double bin;
};

class BDistr {
public:
    BDistr ();
    void Set (ifstream *);
    void Set (double, double);
    void print_param (ostream *);
    double a();
    double b();
private:
    double ain;
    double bin;
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
class NeutronStar {
public:
    parametrs_B * paramet_B;
    double tau; // Время рождения
    bool visible, massive;
    double sparks [2][30];
    double i_incl, x_axis, y_axis, z_axis;
    NeutronStar  (double, OBStar *, MFD *, LM *, GD *, GD *);
    NeutronStar  (double, MFD *, LM *, GD *, GD *, double, double, double, double, double, double, double, double);
    //~neutron_star ();
    double get_M          ();
    double get_R          ();
    double get_P   	(double/*, parametrs_B **/);
    double get_dot_P(double/*, parametrs_B **/);
    //double get_L    (double);
    double get_init_P ();
    double get_init_B ();
    double get_position_x ();
    double get_position_y ();
    double get_position_z ();
    double get_velocity_x ();
    double get_velocity_y ();
    double get_velocity_z ();
    double get_gl (double, SpecialStar *);
    double get_gb (double, SpecialStar *);
    double get_dist_to_sun(double, SpecialStar *);
    void   move_to  (double);
    double get_B    (double/*, parametrs_B **/);
    double get_incl (double);
    //double get_Lum  (double);
    //double get_nu_0 (double);	// Циклическая частота вращения пульсара sec^-1
    //double get_nu_1 (double);	// Первая производная циклической частоты вращения пульсара sec^-2
    //double get_nu_2 (double);	// Вторая производная циклической частоты вращения пульсара sec^-3
    void show_pos_sparks (void);
    void show_pulse_profile (double t, SpecialStar *, parametrs_lum *);
    bool is_pulsar_alive   (double);
    bool is_this_ns();		// Проверка, какая масса у родившейся нейтронной звезды
    double is_pulsar_visible           (double, SpecialStar * , TMap *);
    double get_angle_btw_l_of_s_sparks (double, SpecialStar *, double *, parametrs_lum *);
    double get_DM (double, SpecialStar *, float *, float *, float *);
    //double T_sky  (double, double);
    //double S_min  (double, double, float, double, double, double, float);
private:
MFD * base_mfd;
LM  * base_lm;
double M, R, x, y, z, v_x, v_y, v_z, B, P, L;
};

#endif
