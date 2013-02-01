#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// -----------------------------------------------------------------------
// Функция для реализации численного интегрирования систем ОДУ методом
// Рунге-Кутта-Фельберга 4-5 порядка. Реализована в наиболее обобщенном виде.
// n - количество уравнений
// input[n] - массив начальных значений, результат записывается туда же
// void (*f) (int, double *) - вектор, правая часть уравнений - указатель на эту функцию
// Как показывает опыт эта функция используется только для расчёта небесной механики
// поэтому функция была модифицирована, а именно как минимум 1000 шагов интегрирования
// есть всегда, + 2 шага по 100 000 лет.
// Точность 0.0005 по координатам
// -----------------------------------------------------------------------

void Runge_Kutta (int n, double T, double * input, void (*f)(int n, double * in_fun)) {

    //ofstream out1 ("otl_positionsRKF.txt");


    double result [n], h;
    memset(result, 0, sizeof(result));
    double quant_it;
    double frac_quant, step = 1000000;
    frac_quant = modf (T/step, &quant_it);

    if (quant_it < 0) {
        quant_it *= -1;
        h = -step;
    } else {
        h = step;
    }

    double k_1[n], k_2[n], k_3[n], k_4[n], k_5[n], k_6[n];
    double input1[n], input2[n];
    double vect_for_sent[n];
    double eps;
    double s;
    double t = 0;
    int vn_n=0;

    //cout<<T<<endl;

    for (int i = 0; i < 2*quant_it+1000; i++) 				{
        for (int k = 0; k < n; k++) {
            vect_for_sent[k] = input[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_1[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] + 0.25 *  k_1[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_2[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] + 0.09375 * k_1[k] + 0.28125 * k_2[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_3[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] +  0.879380974 * k_1[k] - 3.27719618 * k_2[k] + 3.32089213 * k_3[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_4[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] +  2.03240741 * k_1[k] - 8 * k_2[k] + 7.17348928 * k_3[k] - 0.205896686 * k_4[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_5[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] - 0.296296296 * k_1[k] + 2 * k_2[k] - 1.38167641 * k_3[k] + 0.45297271 * k_4[k] - 0.275 * k_5[k];
        }

        f(n, &vect_for_sent[0]);
        eps=0;

        for (int k = 0; k < n; k++) 	{
            k_6[k] = h*vect_for_sent[k];
            input1[k] = input[k] + (0.115740741*k_1[k] + 0.548927875*k_3[k] + 0.535331384*k_4[k] - 0.2*k_5[k]);
            input2[k] = input[k] + (0.118518519*k_1[k] + 0.518986355*k_3[k] + 0.50613149 *k_4[k] - 0.18*k_5[k] + 0.036363636*k_6[k]);
            eps += abs(input1[k]-input2[k]);
            //				input[k] += 1./6. * h * (k_1[k] + 2*k_2[k] + 2*k_3[k] + k_4[k]);
        }

        //cout<<"eps - "<<eps<<"\t h - "<<h<<"\t t - "<<t+h<<endl;
        if ((eps<0.0005 && eps>5e-5) || (eps<5e-5 && vn_n>2))	{
            for (int k=0; k < n; k++)	{
                //						out1<<input[k]<<"\t";
                input[k] = input2[k];
            }

            t+=h;

            if (eps>5e-5) {
                vn_n=0;
            }

            //					out1<<endl;
        } else	{
            s = 0.84*pow(0.0005/eps, 0.25);
            //cout<<"s - "<<s<<endl;
            vn_n++;
            h = s*h;
        }

        //				t+=h;
        if (abs(t-T)<h)	{
            //					cout<<"Number of steps - "<<i<<endl;
            break;
        }
    }

    h=T-t;

    quant_it = 1;
    //h = step * frac_quant;
    //if (h<0)
    //h*=sqrt(2.0);


    for (int i = 0; i < quant_it; i++) 				{
        for (int k = 0; k < n; k++) {
            vect_for_sent[k] = input[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_1[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] + 0.25 *  k_1[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_2[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] + 0.09375 * k_1[k] + 0.28125 * k_2[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_3[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] +  0.879380974 * k_1[k] - 3.27719618 * k_2[k] + 3.32089213 * k_3[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_4[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] +  2.03240741 * k_1[k] - 8 * k_2[k] + 7.17348928 * k_3[k] - 0.205896686 * k_4[k];
        }

        f(n, &vect_for_sent[0]);

        for (int k = 0; k < n; k++) 	{
            k_5[k] = h*vect_for_sent[k];
            vect_for_sent[k] = input[k] - 0.296296296 * k_1[k] + 2 * k_2[k] - 1.38167641 * k_3[k] + 0.45297271 * k_4[k] - 0.275 * k_5[k];
        }

        f(n, &vect_for_sent[0]);
        eps=0;

        for (int k = 0; k < n; k++) 	{
            k_6[k] = h*vect_for_sent[k];
            input1[k] = input[k] + (0.115740741*k_1[k] + 0.548927875*k_3[k] + 0.535331384*k_4[k] - 0.2*k_5[k]);
            input2[k] = input[k] + (0.118518519*k_1[k] + 0.518986355*k_3[k] + 0.50613149 *k_4[k] - 0.18*k_5[k] + 0.036363636*k_6[k]);
            eps += abs(input1[k]-input2[k]);
            //				input[k] += 1./6. * h * (k_1[k] + 2*k_2[k] + 2*k_3[k] + k_4[k]);
        }

        //cout<<"eps - "<<eps<<"\t h - "<<h<<"\t t - "<<t+h<<endl;
        //				if (eps<0.0005 && eps>1e-5)			{
        for (int k=0; k < n; k++)	{
            //						out1<<input[k]<<"\t";
            input[k] = input2[k];
        }

        //					out1<<endl;
        //				}
        //				else	{
        //					s = 0.84*pow(0.001/eps, 0.25);
        //cout<<"s - "<<s<<endl;
        //					h = s*h;	}

        t+=h;
        //				if (t<T)
        //					break;
    }





    //cout<<t<<endl;
}

