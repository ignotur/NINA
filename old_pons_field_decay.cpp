#include <cmath>
#include "stars.h"
#include <iostream>
#include <fstream>

using namespace std;

double min (double a, double b) {
	if (a >= b)
		return b;
	else
		return a;
}

//-----------------------------------------------------------------
double MFDOldPons::get_P (double t, double B, double i_incl, double P) {
    double P_res, I, tmp, beta, k1, k2, k3, fractpart, intpart;
    double h_yr, h_sec, sum, alpha, B_min, epart, tau_hall_part;
    double B_0, diff_p, diff_p1, diff_p2, diff_p3, diff_p4;
    double B_diff, diff_t, alpha1, beta1, diff_b, extra_diff_b;
    double lin_approx_b;
    int pos_b, pos_t, int_t;    

    pos_b = 1 + (int) (B / 2.0e11);

    if (pos_b == 0)
	pos_b = 1;
    else if (pos_b > 74)
        pos_b = 74;

    pos_t = (int) t / 1000. - 1;
    int_t = 0; 

    if (pos_t <= 0)
	pos_t = 0;
    else if (pos_t > 10000 - 6) 	{
        int_t = pos_t - 10000;
        pos_t = 10000 - 6;
    }

//	for ( int i=9000; i < 10000; i++)
//		cout << pointer_b[76*i]<<endl; 

    diff_b = B - pointer_b[pos_b];
    extra_diff_b = 0;
    if (pos_b >= 74)	{
        extra_diff_b = B / 1.5e13;
	diff_b = 0;
    }
    diff_t = t - pointer_b[76*pos_t];
    if (pos_t > 10000 - 6)
         diff_t = 0;
   

    diff_p1 = pointer_delta[75 * pos_t + pos_b];
    diff_p2 = pointer_delta[75 * pos_t + pos_b + 1];
    diff_p3 = pointer_delta[75 * (pos_t + 1) + pos_b];
    diff_p4 = pointer_delta[75 * (pos_t + 1) + pos_b + 1];
 
    alpha1 = 1.0 - diff_b/2.0e11; 
    beta1  = 1.0 - diff_t/1000.0;

    diff_p  = beta1 * (alpha1 * diff_p1 + (1 - alpha1) * diff_p2  )  + (1-beta1) * (alpha1 * diff_p3 + (1 - alpha1) * diff_p4);

    if (extra_diff_b == 0)
	    P_res = sqrt(P*P + diff_p + int_t*1000.0 * 3.2e7 * 1.6e-39 * pow(pointer_b[76*pos_t + pos_b], 2.0));
    else
            P_res = sqrt(P*P + diff_p*pow(extra_diff_b, 2.0));     
 
//    cout << "B is "<< B << ", but we used " << pointer_b[pos_b] <<" and pos_b is "<< pos_b<<endl;
//    cout << "Now B is "<< pointer_b[76*pos_t + pos_b] << ", or f(t) = "<< pointer_b[76*pos_t + pos_b] / pointer_b[pos_b]<<endl;
//    cout << "t is "<< t << ", but we used " << pointer_b[76*pos_t] << ", and pos_t is "<< pos_t<< " and int_t is "<< int_t << endl;
//    cout << "diff_p is "<<diff_p<<" and P_res is "<< P_res <<endl;
//    cout << "Additional contribution is "<< int_t*1000.0 * 3.2e7 * 1.6e-39 * pow(pointer_b[76*pos_t + pos_b], 2.0) << endl;
    return P_res;
}

//-------------------------------------------------------------------
double MFDOldPons::get_incl (double t, double B, double i_incl, double P) {
    return i_incl;
}

//-------------------------------------------------------------------
double MFDOldPons::get_dot_P (double t, double B, double i_incl, double P) {
    double res, I, beta = 1.6e-39;
    
    I=1e45;

    res = beta * pow(get_B(t, B, i_incl, P), 2.0) / get_P (t, B, i_incl, P);
    return res;
}
//-------------------------------------------------------------------




//------------------------------------------------------------------

double MFDOldPons::get_B (double t, double B, double i_incl, double P) {
    double res_B, tau_hall_part, B_min, epart;
    int pos_b, pos_t, int_t;    
    double diff_p1, diff_p2, diff_p3, diff_p4, extra_diff_b;
    double diff_p, diff_t, alpha1, beta1, diff_b;

    pos_b = 1 + (int) (B / 2.0e11);

    if (pos_b == 0)
	pos_b = 1;
    else if (pos_b > 74) 
        pos_b = 74;

    pos_t = (int) t / 1000. - 1;
 
    if (pos_t <= 0)
	pos_t = 0;
    else if (pos_t > 10000 - 6) 	{
        int_t = pos_t - 10000;
        pos_t = 10000 - 6;
    }

    diff_b = B - pointer_b[pos_b];
    extra_diff_b = 0;
    if (pos_b >= 74) {
        extra_diff_b = B/1.5e13;
	diff_b = 0;
    }
    diff_t = t - pointer_b[76*pos_t];
//cout << diff_t<<endl;
    if (diff_t > 1000)
         diff_t = 0;
// cout << diff_t<<endl;
  

    diff_p1 = pointer_b[76 * pos_t + pos_b];
    diff_p2 = pointer_b[76 * pos_t + pos_b + 1];
    diff_p3 = pointer_b[76 * (pos_t + 1) + pos_b];
    diff_p4 = pointer_b[76 * (pos_t + 1) + pos_b + 1];
 
    alpha1 = 1.0 - diff_b/2.0e11; 
    beta1  = 1.0 - diff_t/1000.0;

    if (extra_diff_b == 0)
	    diff_p  = beta1 * (alpha1 * diff_p1 + (1 - alpha1) * diff_p2  )  + (1-beta1) * (alpha1 * diff_p3 + (1 - alpha1) * diff_p4);
    else
            diff_p  = extra_diff_b * (beta1 * diff_p1 + (1.0-beta1)*diff_p3);

    res_B = diff_p;
//    res_B = pointer_b[76*pos_t + pos_b];


    return res_B;
}
//------------------------------------------------------------------

void MFDOldPons::print_description (ostream * out) {
    *out<<"#// Model of magnetic field decay is the same as in article  //"<<endl;
    *out<<"#// by Chashkina & Popov (2012) with modification in a       //"<<endl;
    *out<<"#// that first - Hall decay and Coper pair decay, then       //"<<endl;
    *out<<"#// decay because of impurities in crust.                    //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

void MFDOldPons::print_parameters   (ostream * out) {
    *out<<"#//           Parameters of magnetic field decay             //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
    *out<<"#// We read decay model from decay_curve.txt file            //"<<endl;
    *out<<"#//----------------------------------------------------------//"<<endl;
}

MFDOldPons::MFDOldPons (vector <double> * values) {

ifstream in ("decay_curve.txt");

//double  data_per [10000][75];
//double  data_b   [10000][76];

double * data_b;

double * s;

s      = new double [10000*75];
data_b = new double [10000*76];

for (int j = 0; j < 10000; j++)			{
	for (int i = 0; i < 76; i++) {
		if (i==0) 
			in >> data_b   [j*76];
		else {
			in >> data_b   [j*76 + i];
			in >> s [j*75 +  i - 1];
		}
	}
}

//cout << data_b[0][0] << endl;
//cout << data_b[1][0] << endl;
//cout << data_b[2][0] << endl;
//cout << data_b[3][0] << endl;


//cout << s[0]    << endl;
//cout << s[75]   << endl;
//cout << s[2*75] << endl;
//cout << s[3*75] << endl;


//pointer_b     = &data_b[0][0];
pointer_b = data_b;

//for ( int i = 7000; i < 10000; i++ )
//	cout << pointer_b[76*i] << endl;



pointer_delta = s; 
}

