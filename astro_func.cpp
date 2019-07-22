#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include "stars.h"
#include <cstdlib>

using namespace std;

//--------------------------------------------------------------//
// Функция, которая возвращает температуру неба в направлении
// детектируемого пульсара на частоте 1.4 МГц.
// Основана на работе Dinnat et al., 2010
// Аргументы функции - галактическая долгота и широта пульсара
//--------------------------------------------------------------//

bool isNaN(double x) { 
  return x != x;
}

struct decl {
    double quant;
};

double T_sky (double l, double b, TMap * T_copy) {

    double alphaf, deltaf, res;
    double alpha_int, alpha_frac, delta_int, delta_frac;
    double N_deg;
    double A_g [3][3], ecv[3], gal[3];
    decl Tb;

    //ifstream in_dec ("declination.bin", ios::binary);
    //ifstream in_ras ("right_ascension.bin", ios::binary);
    //ifstream in_Tb  ("TbGal_tot_CasA1pix.bin", ios::binary);

    A_g [0][0] = -0.0548755601367195;
    A_g [0][1] = +0.49410942801324;
    A_g [0][2] = -0.86766614895829;
    A_g [1][0] = -0.87343709025327;
    A_g [1][1] = -0.4448296298016944;
    A_g [1][2] = -0.19807637370567;
    A_g [2][0] = -0.48383501554722;
    A_g [2][1] = +0.74698224450044;
    A_g [2][2] = +0.4559837761713720;

    //double l = 20./180.*pi,  b = 50./180.*pi;

    //cout<<l<<"\t"<<b<<endl;

    gal[0] = cos(b)*cos(l);
    gal[1] = cos(b)*sin(l);
    gal[2] = sin(b);

    memset(ecv, 0, sizeof(ecv));

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            ecv[i] += A_g[i][j]*gal[j];
        }

    deltaf = asin(ecv[2]);
    alphaf = atan2(ecv[1], ecv[0]);

    deltaf*=180./pi;
    alphaf*=180./pi;

    if (alphaf<0) {
        alphaf+=360;
    }

    //cout<<alphaf<<"\t"<<deltaf<<endl;

    alphaf*=4.;
    deltaf+=90.;
    deltaf*=4.;

    alpha_frac = modf (alphaf, &alpha_int);
    delta_frac = modf (deltaf, &delta_int);



    N_deg = 721*alpha_int;
    N_deg += delta_int;

    //	in_dec.seekg(sizeof(tmp)*N_deg, ios_base::beg);
    //	in_dec.read((char*) &delta, sizeof(tmp));
    //	in_ras.seekg(sizeof(tmp)*N_deg, ios_base::beg);
    //	in_ras.read((char*) &alpha, sizeof(tmp));
    //	in_Tb.seekg(sizeof(Tb)*N_deg, ios_base::beg);
    //	in_Tb.read((char*) &Tb, sizeof(Tb));

    res = T_copy->get_Tb(N_deg);

    //cout<<l<<"\t"<<b<<"\t"<<res<<endl;

    //res = Tb.quant;
    return res;
}

double S_min (double l, double b, float sm, double dist, double w, double P, float DM, TMap * T_copy) {
    double res, DM_0_parkes, DM_0_swinburne, delta_beam, tau_scatt;
    double W_l_parkes, W_l_swinburne, S_min_Parkes, S_min_Swinburne;

    double N_ch = 96, t_sampl_parkes = 250e-6, t_sampl_swinburne = 125e-6, nu = 1.4e3, delta_nu = 288e6;  // nu and delta_nu should be measured in MHz in principle
    double tau_sampl_parkes = 1.5*t_sampl_parkes, tau_sampl_swinburne = 1.5*t_sampl_swinburne;
    double beta = 1.5, sigma = 8, T_rec = 21, Tb_sky, G = 0.64, N_p = 2, t_int_parkes = 2100.0, t_int_swinburne = 265.0; // G = 0.64 K/Jy

//    tau_scatt = 1000.*pow(sm/292., 1.2)*dist*pow(nu, -4.4);

    //if (dist<20)
    //cout<<"DM - "<<DM<<", dist - "<<dist<<", DM/dist - "<<DM/dist<<endl;

    //------------------------------------------------------//
    // Lorimer et al. ArXiv:0607640
    //------------------------------------------------------//

//    tau_scatt = 0.154*log10(DM)+1.07*pow(log10(DM), 2.) - 7.;
  //  tau_scatt = pow(10., tau_scatt)/1.e3;

    //------------------------------------------------------//

    //------------------------------------------------------//
    //  Bhat, Cordes, Camilo et al. (2004)
    //------------------------------------------------------//

    tau_scatt = pow(10.0, (-6.46 + 0.154 * log10(DM) + 1.07 * pow(log10(DM), 2))) * pow((nu/1e3), -3.86)/1000.0;



    delta_beam = exp(-pow(rand() * sqrt(2.0*log(2.0)) /rand_high_board, 2));
   // delta_beam = 0.5;
    Tb_sky = T_sky (l,b, T_copy);

    DM_0_parkes = N_ch*t_sampl_parkes*pow(nu,3)/8299./(delta_nu/1.e6);
    DM_0_swinburne = N_ch*t_sampl_swinburne*pow(nu,3)/8299./(delta_nu/1.e6);
    W_l_parkes = sqrt(w*w + tau_sampl_parkes*tau_sampl_parkes + pow(t_sampl_parkes*DM/DM_0_parkes, 2) + tau_scatt*tau_scatt);
    W_l_swinburne = sqrt(w*w + tau_sampl_swinburne*tau_sampl_swinburne + pow(t_sampl_swinburne*DM/DM_0_swinburne, 2) + tau_scatt*tau_scatt);

    if (P > W_l_parkes) {
	    S_min_Parkes = delta_beam * beta*sigma*(T_rec + Tb_sky)/G/sqrt(N_p*delta_nu*t_int_parkes)*sqrt(W_l_parkes/(P-W_l_parkes));
	    S_min_Swinburne = delta_beam * beta*sigma*(T_rec + Tb_sky)/G/sqrt(N_p*delta_nu*t_int_swinburne)*sqrt(W_l_swinburne/(P-W_l_swinburne));
    }
    else   {
	    S_min_Parkes = 1e9;
	    S_min_Swinburne = 1e9;
    }


    if (abs(b*180./3.1415926) >= 5.) {
        res = S_min_Swinburne;
    } else {
        res = S_min_Parkes;
    }


//    cout << "S_min_Parkes -- "<<S_min_Parkes << ", S_min_Swinburne -- "<<S_min_Swinburne << endl;
    if (isNaN(S_min_Parkes))	{
	cout << "Debugging -- " << P - W_l_parkes << "\t" << P << "\t" << w << "\t" <<sqrt(tau_scatt*tau_scatt)<<endl; 
	cout << "DM -- "<<DM<<endl;
    }
    //res = min (S_min_Parkes, S_min_Swinburne);

    return res;
}
