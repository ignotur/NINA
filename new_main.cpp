#include <iostream>
#include "stars.h"
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstring>

using namespace std;

int main (int argc, char * argv[]) {

    ifstream input     ("input.txt");

    parametrs_B   param_B   (&input);
    parametrs_lum param_lum (&input);

    double T = -5.e8; //-350e6;              // начало
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

    PDistr p_distr;
    BDistr b_distr;	

    if (argc>2)				{
        in_p.open (argv[1]);
        in_b.open (argv[2]);
    } else										{
        cout<<"You did not mention files for P and B distributions, so standards files are used"<<endl;
        in_p.open ("P_init.txt");
        in_b.open ("B_init.txt");
		
	if ((!in_p) || (!in_b))						{
	    cout<<"Files with initial P&B are absent."<<endl;
	} else 							{		
	    p_distr.Set(&in_p);
	    b_distr.Set(&in_b);
	}		

    }



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
    ofstream out_ances ("result_ancest_pos.txt");
    ofstream out_about ("output.txt");
    ofstream out_short ("output_short.txt");

    ofstream out_exper ("result_experiment.txt");

    ofstream out_str ("otl_v_x.txt");

    param_B.print_short    (&out_short);
    param_lum.print_short  (&out_short);



    SpecialStar sun;
    SpecialStar sun_nowaday;
    double now = 0, shift;
    double P, dot_P, x, y, z, B;
    double dist_to_sun, lumin;
    double x_ances, y_ances, z_ances;
    int counter = 0;
    int rand_shift;                 // случайное время рождения внутри тысячелетия
    int n_magnet = 0;
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
    p_distr.print_param        (&out_about);
    b_distr.print_param        (&out_about);
    out_about<<"//----------------------------------------------------------//"<<endl;
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    out_about<<"// Запушено "<<asctime(timeinfo);


    OBStar      * ancester;
    NeutronStar * descendant;

    for (int i = 0; i < number_millenium; i++)		{
        sun.move_to(T);

        for (int j = 0; j < number_stars; j++)	{
            rand_shift = rand()%500;
            ancester = new OBStar (T + rand_shift, &sun);

            if (!(i%10000) && j==0) {
                cout<<T<<"\t"<<i<<endl;
            }

            shift  = ancester->get_time_on_MS();
            shift += ancester->get_t_He      ();

            if (shift + T + rand_shift < 0)		{
                ancester->move_to (shift); // На сколько нужно сдвинуть, действительно мы не знаем времени рождения звезды
                descendant = new NeutronStar (T + shift + rand_shift, ancester, &param_B, &p_distr, &b_distr);
                P     = descendant->get_P(now/*,&param_B*/);
                out_str << descendant->get_velocity_x()<<endl;

                //			cout<<P<<endl;
                if (descendant->is_pulsar_alive(now) &&  descendant->is_this_ns() && param_lum.is_beam_on(P))	{
                    //cout<<"TTT"<<endl;
                    //cout<<descendant->get_B(now)<<"\t"<<descendant->get_dot_P(now)<<endl;
                    descendant->move_to(now);
                    lumin = descendant->is_pulsar_visible(now, &sun_nowaday, &T_copy, &param_lum);

                    if (lumin)				{
                        counter++;
                        //						if (descendant->is_this_ns())	{
                        //						x_ances  = ancester->get_position_x();
                        //						y_ances  = ancester->get_position_y();
                        //						out_ances<<sqrt(pow(x_ances, 2) + pow(y_ances, 2))<<endl;	}
                        x     = descendant->get_position_x();
                        y     = descendant->get_position_y();
                        z     = descendant->get_position_z();
                        P     = descendant->get_P(now/*,&param_B*/);
                        dot_P = descendant->get_dot_P(now/*,&param_B*/);
                        B     = descendant->get_B (now);
                        dist_to_sun = descendant->get_dist_to_sun(now, &sun_nowaday);

                        if (dist_to_sun < 10) {
                            out_p  << P<<"\t"<< dot_P<<endl;
                        }

                        out_pos<< x<<"\t" << y <<"\t"<< z <<endl;
                        out_lum<<lumin*1000*pow(dist_to_sun, 2)<<endl;

                        //						out_exper<<P<<"\t"<<P/2./dot_P/365.24/24/3600.<<"\t"<<-(T + shift + rand_shift)<<endl;
                        out_exper<<P*sqrt(1 + (T + shift + rand_shift)/(P/2./dot_P/365.24/24/3600.))<<endl;

                        if (B > 1e14) {
                            n_magnet++;
                        }
                    }
                }

                delete descendant;
            }

            delete ancester;
        }

        T += 1000;
    }

    time (&rawtime);
    timeinfo = localtime(&rawtime);
    out_about<<"// Закончено "<<asctime(timeinfo);
    out_about<<"// За время работы пульсаров замечено  "<<counter<<endl;
    out_about<<"// За время работы магнетаров замечено "<<n_magnet<<endl;

    return 0;
}
