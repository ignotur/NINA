#include <cstdlib>
#include <time.h>
#include <cmath>
#include <iostream>
#include "stars.h"
#include <algorithm>
#include <cstring>

using namespace std;

//-------------------------------------------------------------------------------------//
// neutron_star --- класс, который на месте массивной звезды создаёт нейтронную звезду,
// а также позволяет прослеживать её дальнейшую эволюцию. Конструктору класса передаётся
// звезда предшественник этой нейтронной звезды. Извлекаеться положение в пространстве
// Галактике, магнитный поток, масса ядра прародителя
// Начальная скорость вращения взята из статьи Бисноватого-Когана, 2008
// Масса по статье Hurley, 2000
//
//  Автор: Igoshev Andrey
//  Научный руководитель: Холтыгин А.Ф.
//  e-mail: igoshev-andrei@rambler.ru
//  Написание начато: 13.09.2010
//-------------------------------------------------------------------------------------//

//------------------------------------------------------------------//
// Декларация функций
double rho_hartman (double);
double norm_distr  (void);
double expon_vel   (double);
void   Runge_Kutta (int, double, double *, void (*f)(int, double *));
void   diff_equi   (int, double *);
//------------------------------------------------------------------//

neutron_star::neutron_star (double T, star_OB  * proteg, parametrs_B * param, P_distr *p_init, B_distr *b_init) {
    double v, prover, sigma, phi, psi, f, chance_1, chance_2, chance_3;
    bool is_position_set = false;

    M = 1.17 + 0.09 * proteg->get_M_c_SN();

    // Определение радиуса нейтронной звезды -------------------
    // по статье J.M. Latimer, 2000
    R = 3.04*G_cgs*M*M_sol/pow(light_velocity, 2);
    R = R / sqrt(1-2*G_cgs*M*M_sol/(R*pow(light_velocity,2)));
    //-----------------------------------------------------------

    // Начальный период пульсара распределён с математическим
    // ожиданием 6 миллисекунд и дисперсией 2 миллисекунды
    // оценка по статье Бисноватый-Коган, 2008 АЖ 12
    //P = 6e-3 + 2e-3*norm_distr();
    //<--------Только в целях отладки и тестирования
    //chance_1  = rand()/rand_high_board;
    //chance_1 *= 3;

    //	if (chance_1>1.)
    //		P = 0.02 + 0.02*norm_distr();
    //	else
    //		P = 0.2  + 0.05*norm_distr();
    //P = 0.3 + 0.15*norm_distr();
    //P = 0.1/* + 0.1*norm_distr()*/;
    P=p_init->a()+p_init->b()*norm_distr();
    //cout<<P<<endl;
    //P = 0.15 + 0.08*norm_distr();
    //P = 0.08 + 0.06*norm_distr();
    //P = 0.26;
    //P=0.008;
    //------------------------------------------------------------

    // Оценка кинетической энергии теряемой пульсаром при сбросе
    // оболочки за счёт эффектов магниторотационного взрыва
    L =  proteg->get_L();
    L -= 4./5.*M*M_sol*pow(R,2)*pi/P;
    L /= 0.2;
    // В статье Бисноватого-Когана процесс продолжается
    // около 0,2 секунд, всё изменение момента за это время
    // переходит в энергию.
    //------------------------------------------------------------

    // Координаты звезды наследуются от её предшественника
    x = proteg->get_position_x();
    y = proteg->get_position_y();
    z = proteg->get_position_z();

    //-------------------------------------------------------------
    // Скорость звезды меняеться на начальный толчок (birthkick)
    // распределённый случайным образом с плотностью вероятности
    // распределённой по нормальному закону
    double is_velocity_set;
    //v_z = 150*1e5*lcm/lsec*norm_distr() + proteg->get_velocity_x();
    //v_x = 150*norm_distr()*1e5*lcm/lsec + proteg->get_velocity_y();
    //v_y = 150*norm_distr()*1e5*lcm/lsec + proteg->get_velocity_z();
    //Повторим Faucher для начала
    /*is_velocity_set = false;
    	do {
    		chance_1 = rand () / rand_high_board;
    		chance_1 *= 1200;
    		chance_1 += 600;
    		chance_2 = rand () / rand_high_board;

    		if (expon_vel(chance_1) >= chance_2)
    			is_velocity_set = true;

    	} while (!(is_position_set));

    cout<<"first"<<endl;

    v_x = chance_1*1e5*lcm/lsec;
    is_velocity_set = false;
    	do {
    		chance_1 = rand () / rand_high_board;
    		chance_1 *= 1200;
    		chance_1 += 600;
    		chance_2 = rand () / rand_high_board;

    		if (expon_vel(chance_1) >= chance_2)
    			is_velocity_set = true;

    	} while (!(is_position_set));

    cout<<"second"<<endl;

    v_y = chance_1*1e5*lcm/lsec;
    is_velocity_set = false;
    	do {
    		chance_1 = rand () / rand_high_board;
    		chance_1 *= 1200;
    		chance_1 += 600;
    		chance_2 = rand () / rand_high_board;

    		if (expon_vel(chance_1) >= chance_2)
    			is_velocity_set = true;

    	} while (!(is_position_set));
    v_z = chance_1*1e5*lcm/lsec;
    */
    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_x = 1e5*lcm/lsec*expon_vel(chance_1);
    } else {
        v_x = -1e5*lcm/lsec*expon_vel(chance_1);
    }

    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_y = 1e5*lcm/lsec*expon_vel(chance_1);
    } else {
        v_y = -1e5*lcm/lsec*expon_vel(chance_1);
    }

    chance_1 = rand () / rand_high_board;
    chance_1/= 360.;
    chance_2 = rand () / rand_high_board;

    if (chance_2 > 0.5) {
        v_z = 1e5*lcm/lsec*expon_vel(chance_1);
    } else {
        v_z = -1e5*lcm/lsec*expon_vel(chance_1);
    }

    //ofstream output_str("otl_v_x.txt");
    //output_str<<expon_vel(chance_1)<<endl;

    //-------------------------------------------------------------

    //-------------------------------------------------------------
    // При вспышке сверхновой уноситься магнитное поле пропорциональное
    // унесённой массе. Физическое обоснование - вмороженность
    // силовых линий поля в плазму
    //B = proteg->get_Flux();
    //B = M/proteg->get_mass()*B;
    //B = B/4./pi/pow(R, 2);
    //<-----------Только для тестирования!

    //chance_1  = rand()/rand_high_board;
    //chance_1 *= 5;

    //	if (chance_1>1.)
    //		B = 12.20+0.35*norm_distr();
    //	else
    //		B = 13.30+0.35*norm_distr();
    //B = 13.25+0.6*norm_distr();
    //B=12.65+0.55*norm_distr();
    //B = 13.2+0.7*norm_distr();
    B=b_init->a()+b_init->b()*norm_distr();
    B=pow(10,B);
    //cout<<B<<endl;
    //B = 37e12;
    //B=1e11;
    //--------------------------------------------------------------

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
    //cout<<"intcl "<<i_incl*180/pi<<endl;

    //-----------------------------------------------------------------
    // Устанавливаем текущее время временем рождения нейтронной звезды
    tau = T;
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    // Не слишком ли массивна эта звезда, чтобы быть нейтронной?
    if (proteg->get_c_BAGB()<9.5084) {
        massive = true;
    } else {
        massive = false;
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


    //parametrs paramet_B (param);
    paramet_B = param;
    //-----------------------------------------------------------------
}

double neutron_star::get_M() {
    return M;
}

double neutron_star::get_R() {
    return R/1e5;
}
/*
//-----------------------------------------------------------------
// Период вращения нейтронной звезды. Формулы для периода с учётом
// убывания магнитного поля получены в основной работе методом
// интегрирования дифференциального уравнения для электромагнитного
// дипольного излучения и токовых потерь

double neutron_star::get_P(double t) {
double P_res, I, tmp;
I = 2./5. * M*M_sol*pow(R,2);
t = t - tau;
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

//-------------------------------------------------------------------
// Первая производная периода. Формулы взяты из исходного дифферен-
// циального уравнения

double neutron_star::get_dot_P (double t) {
double res, I;
t = t - tau;
I = 2./5. * M*M_sol*pow(R,2);
double tau_mu = 3.*pow(light_velocity,3)*I/2./pow(B,2)/pow(R,6);
	if (t>1e6)
	res = pow(e, -2.*1e6/tau_ohm)/get_P(t+tau)/tau_mu/pow(1.+alpha*(1-pow(e, -1e6/tau_ohm)),2);
	else
	res = pow(e, -2.*t/tau_ohm)/get_P(t+tau)/tau_mu/pow(1.+alpha*(1-pow(e, -t/tau_ohm)),2);
return res;
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// Следующие функции возвращают первую, вторую производные и саму
// частоту вращения нейтронной звезды.
// Вторая производная частоты вращения взята из дифференциального
// уравнения.

double neutron_star::get_nu_0 (double t) {
double res;
res = 1./get_P (t);
return res;
}

double neutron_star::get_nu_1 (double t) {
double res;
res = - 1.*pow(get_nu_0(t), 2) * get_dot_P (t);
return res;
}

double neutron_star::get_nu_2 (double t) {
double res;
	if (t-tau<1e6)
	res = - 2.*get_nu_1(t) / tau_ohm / lsec + pow(get_nu_1(t), 2) / get_nu_0(t) - alpha * pow (e, -(t-tau)/tau_ohm) * get_nu_1(t) / tau_ohm / lsec + 2. * get_nu_1(t) / get_nu_0(t);
	else
	res = pow(get_nu_1(t),2) / get_nu_0(t) + 2.*get_nu_1(t) / get_nu_0(t);
return res;
}
*/
//---------------------------------------------------------------------

/*
//---------------------------------------------------------------------
// Модель светимость - функция описываюзия зависимость от возраста
// пульсара

double neutron_star::get_Lum (double t)  {
double res;

t = t - tau;
if (t>1000)			{
res = -0.646*log10(t) + 7.524;
res = pow(10, res);		}
else
res = 1000;
return res;
}

//-----------------------------------------------------------------------
*/
//-----------------------------------------------------------------------
// Не погас ли ещё пульсар? Алгоритм проверки как у Faucher

bool neutron_star::is_pulsar_alive(double t) {

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

double neutron_star::get_position_x() {
    return x;
}

double neutron_star::get_position_y() {
    return y;
}

double neutron_star::get_position_z() {
    return z;
}

void neutron_star::move_to(double T) {
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

double neutron_star::get_velocity_x (void) {
    return v_x/lcm*lsec/1e5;
}

double neutron_star::get_velocity_y (void) {
    return v_y/lcm*lsec/1e5;
}

double neutron_star::get_velocity_z (void) {
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
bool neutron_star::is_this_ns (void) {
    return massive;
}

double neutron_star::get_dist_to_sun(double t, special_star * sun)	{
    double res;
    double x_sun, y_sun, z_sun;

    sun->move_to(t);

    x_sun = sun->get_position_x();
    y_sun = sun->get_position_y();
    z_sun = sun->get_position_z();

    res = sqrt(pow(x - x_sun, 2) + pow(y - y_sun, 2) + pow(z - z_sun, 2));

    return res;
}

//---------------------------------------------------------------//
// Получение меры дисперсии для данного наблюения пульсара,
// используется код NE2001
//---------------------------------------------------------------//

extern "C" {
    void dmdsm_ (float *l, float *b, int *ndir, float *dmpsr, float *dist, char *limit, float *sm, float *smtau, float *smtheta, float *smiso);
}


double neutron_star::get_DM (double t, special_star * sun, float *l, float *b, float *sm) {
    double x_sun, y_sun, z_sun;
    float dist, dmpsr;
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

    int ndir = -5;

    float /*sm, */ smtau, smtheta, smiso;

    char limit;
    limit = ' ';

    //cout<<*l<<"\t"<<*b<<"\t"<<dist<<endl;

    dmdsm_ (l, b, &ndir, &dmpsr, &dist, &limit, sm, &smtau, &smtheta, &smiso);

    return dmpsr;
}
