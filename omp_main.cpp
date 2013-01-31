#include <iostream>
#include "stars.h"
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

int main (int argc, char *argv[]) {

    ifstream input     ("input.txt");

    parametrs_B   param_B   (&input);
    parametrs_lum param_lum (&input);

    double T = -5e8;              // начало
    int number_stars = 7;		// темп звездообразования (звёзд в тысячалетие)
    int number_millenium = -T/1e3;  // количество тысячалетий

    cout<<"//----------------------------------------------------------//"<<endl;
    cout<<"// Пульсарная популяция. Версия 0.80см соглашение о версиях //"<<endl;
    cout<<"// Автор: Игошев Андрей, научный руководитель: А.Ф. Холтыгин//"<<endl;
    cout<<"// e-mail: igoshev-andrei@rambler.ru СПбГУ, 2010-2012       //"<<endl;
    cout<<"//----------------------------------------------------------//"<<endl;
    param_B.print_description  (&cout);
    param_B.print_parametrs    (&cout);
    param_lum.print_description(&cout);
    param_lum.print_parametrs  (&cout);

    ifstream in_p, in_b;

    if (argc>2)				{
        in_p.open (argv[1]);
        in_b.open (argv[2]);
    } else										{
        cout<<"You did not mentioned files for P and B distributions, so standards files are used"<<endl;
        in_p.open ("P_init.txt");
        in_b.open ("B_init.txt");
    }

    PDistr p_distr(&in_p);
    BDistr b_distr(&in_b);

    p_distr.print_param        (&cout);
    b_distr.print_param        (&cout);

    cout<<"// Параметры популяциии:                                          "<<endl;
    cout<<"// T_start "<<T<<endl;
    cout<<"// star formation rate "<<number_stars<<endl;
    cout<<"//----------------------------------------------------------//"<<endl;
    cout<<"// Инициализация расчётов.                                  //"<<endl;
    srand(time(0));
    TMap T_copy;

    cout<<"// Инициализация расчётов закончена.                        //"<<endl;
    cout<<"//----------------------------------------------------------//"<<endl;


    ofstream out_p     ("result_p_p_dot.txt");
    ofstream out_pos   ("result_pos.txt");
    ofstream out_lum   ("result_radio_lum.txt");
    ofstream out_about ("output.txt");
    ofstream out_short ("output_short.txt");

    param_B.print_short    (&out_short);
    param_lum.print_short  (&out_short);


    SpecialStar sun;
    SpecialStar sun_nowaday;
    double now = 0, shift;
    double P[number_stars], dot_P[number_stars], x[number_stars], y[number_stars], z[number_stars], B[number_stars];
    double dist_to_sun[number_stars], lumin;
    int counter = 0;
    int rand_shift;                 // случайное время рождения внутри тысячелетия
    int n_magnet = 0;
    bool active[number_stars];
    sun_nowaday.move_to(now);

    out_about<<"//----------------------------------------------------------//"<<endl;
    out_about<<"// Пульсарная популяция. Версия 0.80см соглашение о версиях //"<<endl;
    out_about<<"// Автор: Игошев Андрей, научный руководитель: А.Ф. Холтыгин//"<<endl;
    out_about<<"// e-mail: igoshev-andrei@rambler.ru СПбГУ, 2010-2012       //"<<endl;
    out_about<<"//----------------------------------------------------------//"<<endl;
    out_about<<"// Параметры модели:                                          "<<endl;
    out_about<<"// T_start "<<T<<endl;
    out_about<<"// star formation rate "<<number_stars<<endl;
    param_B.print_description  (&out_about);
    param_B.print_parametrs    (&out_about);
    param_lum.print_description(&out_about);
    param_lum.print_parametrs  (&out_about);
    out_about<<"//----------------------------------------------------------//"<<endl;
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    out_about<<"// Запушено "<<asctime(timeinfo);

    int rank[7];

    OBStar      * ancester;
    NeutronStar * descendant[number_stars];

#ifdef _OPENMP
    omp_set_num_threads(3);
#endif

    for (int i = 0; i < number_millenium; i++)		{
        sun.move_to(T);

        memset (active, false, sizeof(active));

        #pragma omp parallel private(rand_shift, ancester, shift) shared(T, n_magnet, counter, rank, descendant, x, y, z, P, dot_P, B, dist_to_sun, active)
        {
            #pragma omp for

            for (int j = 0; j < number_stars; j++)	{
                rand_shift = rand()%500;
                ancester = new OBStar (T + rand_shift, &sun);

                if (!(i%10000) && j==0) {
                    cout<<T<<"\t"<<i<<endl;
                }

                shift  = ancester->get_time_on_MS();
                shift += ancester->get_t_He      ();


                if (shift + T + rand_shift < 0)		{
                    active[j]=true;
                    ancester->move_to (shift); // На сколько нужно сдвинуть, действительно мы не знаем времени рождения звезды
                    descendant[j] = new NeutronStar (T + shift + rand_shift, ancester, &param_B, &p_distr, &b_distr);

                    if (descendant[j]->is_pulsar_alive(now) && descendant[j]->is_this_ns())	{
                        descendant[j]->move_to(now);
                        //	lumin = descendant->is_pulsar_visible(now, &sun_nowaday, &T_copy, &param_lum);
                        //		if (lumin)				{
                        //			counter++;
                        x[j]     = descendant[j]->get_position_x();
                        y[j]     = descendant[j]->get_position_y();
                        z[j]     = descendant[j]->get_position_z();
                        P[j]     = descendant[j]->get_P(now/*,&param_B*/);
                        dot_P[j] = descendant[j]->get_dot_P(now/*,&param_B*/);
                        B[j]     = descendant[j]->get_B (now);
                        dist_to_sun[j] = descendant[j]->get_dist_to_sun(now, &sun_nowaday);
                        //
                        //cout<<B[j]<<endl;
                        //			rank[j]=omp_get_thread_num()+1;
                        //			cout<<rank<<"\t"<<j<<endl;
                        //			out_p  << P<<"\t"<< dot_P<<endl;
                        //			out_pos<< x<<"\t" << y <<"\t"<< z <<endl;
                        //			out_lum<<lumin*1000*pow(dist_to_sun, 2)<<endl;
                        //}
                        //						if (B > 1e14)
                        //							n_magnet++;		}
                    }

                    //	delete descendant;
                } else	{
                    //				cout<<"first sigh"<<endl;
                    active[j]=false;
                }

                delete ancester;
            }
        }
        T += 1000;

        for (int j = 0; j < number_stars; j++)	{
            if (active[j])
                if (descendant[j]->is_pulsar_alive(now) && descendant[j]->is_this_ns()) {
                    lumin = descendant[j]->is_pulsar_visible(now, &sun_nowaday, &T_copy, &param_lum);

                    if (lumin)	{
                        counter++;
                        out_p  << P[j]<<"\t"<< dot_P[j]<<endl;
                        out_pos<< x[j]<<"\t" << y[j] <<"\t"<< z[j] <<endl;
                        out_lum<<lumin*1000*pow(dist_to_sun[j], 2)<<endl;

                        if (B[j] > 1e14) {
                            n_magnet++;
                        }


                    }
                }

            if (active[j]) {
                delete descendant[j];
            }
        }


    }

    time (&rawtime);
    timeinfo = localtime(&rawtime);
    out_about<<"// Закончено "<<asctime(timeinfo);
    out_about<<"// За время работы пульсаров замечено  "<<counter<<endl;
    out_about<<"// За время работы магнетаров замечено "<<n_magnet<<endl;

    return 0;
}
