#include <cmath>
#include "stars.h"
#include <iostream>
#include <fstream>

using namespace std;

//--------------------------------------------------------------//
// Этот файл содержит конструктор и методы класса
// special_star - выделенная звезда, в первоначальном
// проекте только солнце, для определения параметров
// спиральных рукавов в любой момент времени,
// и в последующем для проведения "наблюдений"
// пульсарной популяции.
// Автор: Igoshev Andrey
// Написание начато: 2.10.2010
//---------------------------------------------------------------//

//---------------------------------------------------------------//
// Декларация функций
void Runge_Kutta (int, double, double *, void (*f)(int, double *));
void diff_equi   (int, double *);
double dphi_dr (double, double, double);
double dphi_dx (double, double, double);
double dphi_dy (double, double, double);
//---------------------------------------------------------------//

//---------------------------------------------------------------//
// Скорость солнца в галактике взята из книги Аллен
SpecialStar::SpecialStar() {
    double r, v;
    x = 0;
    y = 8.5;
    z = 0;
    r = sqrt(x*x+y*y+z*z);

    v_x = 5.23*1e5*lcm/lsec;
    v_y = -10.*1e5*lcm/lsec;
    v_z = 7.17*1e5*lcm/lsec;
    //v_x=0;
    //v_y=0;
    //v_z=0;


    v = sqrt(r*dphi_dr(x,y,0));
    v_x += v;
    t=0;
    use = 0;

    if (std::ifstream("sun_movement.txt") != NULL) {
        out_err.open("sun_movement_add_sun.txt");
    } else {
        out_err.open("sun_movement.txt");
    }
}


double SpecialStar::get_position_x () {
    return x;
}

double SpecialStar::get_position_y () {
    return y;
}

double SpecialStar::get_position_z () {
    return z;
}

double SpecialStar::get_velocity_x () {
    return v_x/lcm*lsec/1e5;
}

double SpecialStar::get_velocity_y () {
    return v_y/lcm*lsec/1e5;
}

double SpecialStar::get_velocity_z () {
    return v_z/lcm*lsec/1e5;
}

//---------------------------------------------------------//
// Метод позволяющий вычислять положение солнца в любой
// момент времени. Время - абсолютно. Солнце рождается
// в современной Галактике, а потом движеться назад, если то
// необходимо.
// Метод составляет диффур движения с начальными параметрами
// и вызывает интегратор.
//---------------------------------------------------------//
void SpecialStar::move_to(double T) {
    double result [6];
    double r, v;

    out_err<<T-t<<endl;

    if (use>5000)	{
        //cout<<"reset"<<endl;
        x = 0;
        y = 8.5;
        z = 0;
        r = sqrt(x*x+y*y+z*z);

        v_x = 5.23*1e5*lcm/lsec;
        v_y = -10.*1e5*lcm/lsec;
        v_z = 7.17*1e5*lcm/lsec;

        v = sqrt(r*dphi_dr(x,y,0));
        v_x += v;
        use = 0;
        t=0;
    }


    if (t != T)	{
        T = T-t;
        result [0] = x;
        result [1] = y;
        result [2] = z;
        result [3] = v_x;
        result [4] = v_y;
        result [5] = v_z;

        //	cout<<"The sun is moving to "<<T<<endl;

        Runge_Kutta (6, T, &result[0], &diff_equi);

        x   = result [0];
        y   = result [1];
        z   = result [2];
        v_x = result [3];
        v_y = result [4];
        v_z = result [5];
        t = t + T;
        use++;
    }

    /*
    	else if (t==0)	{
    cout<<"reset"<<endl;
    		x = 0;
    		y = 8.5;
    		z = 0;
    		r = sqrt(x*x+y*y+z*z);

    		v_x = 5.23*1e5*lcm/lsec;
    		v_y = -10.*1e5*lcm/lsec;
    		v_z = 7.17*1e5*lcm/lsec;

    		v = sqrt(r*dphi_dr(x,y,0));
    		v_x += v;
    	t=0;		}

    */
}

//-----------------------------------------------//
// Функция возвращает позиционный угол солнца
// в галактике. Нужна для определения положения
// спиральных рукавов в любой момент времени
//-----------------------------------------------//
double SpecialStar::get_theta() {
    double alpha;
    alpha = atan2(x, y);

    if (alpha < 0) {
        alpha += 2*pi;
    }

    return alpha;
}
