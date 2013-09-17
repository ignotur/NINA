#include <iostream>
#include "stars.h"
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <unistd.h>

using namespace std;

int main (int argc, char * argv[]) {

print_head(&cout);

int rez=0, type_of_run;
// There are following type of run here:
// 0 - error, 2 - full start, 4 - generate massive stars,
// 5 - use already generated massive stars
bool ext_print =false, file_for_print=false, hidden_file = false;
bool black_hole=false; 
bool exclusion = false, param_file=false;
char * name_file_param, * argum, * file_res, * read_file, * BH_file;
type_of_run = 0;

    while ( (rez = getopt(argc,argv,"hefo:i:g:p:ub:")) != -1){
        switch (rez){
	case 'h': type_of_run = 0   ; break;
	case 'e': ext_print   = true; break;
	case 'f': if (type_of_run == 4 || type_of_run == 5) exclusion=true; else type_of_run = 2; break;
	case 'o': file_for_print = true; file_res = optarg; break;
	case 'g': if (type_of_run == 2 || type_of_run == 5 || black_hole) exclusion=true; else type_of_run = 4; argum     = optarg; break;
	case 'i': if (type_of_run == 2 || type_of_run == 4) exclusion=true; else type_of_run = 5; read_file = optarg; break;
	case 'p': param_file  = true; name_file_param = optarg; break;
	case 'b': if (type_of_run == 4) exclusion=true; else { black_hole = true; BH_file = optarg; } break;
	case 'u': hidden_file = true; break;
	case '?': type_of_run = 0; break;
	default : type_of_run = 4; break;
   };
};

    if (type_of_run == 0)	{
	print_help();
	return 0;
    }

    if (exclusion)	{
	print_exclusion();
	return 1;
    }	 

ifstream in_param;
bool arms;
char star_distr, init_distr_p, init_distr_f, lum_model, decay;
vector <double> param_p, param_f, list_lum, param_decay;
double time_of_run, birthrate;

cout<<"param_file is "<<param_file<<endl;

    if (param_file)	{
	in_param.open (name_file_param);
	input_syntax (&in_param, &star_distr, &arms, &init_distr_p, &param_p, 
                   &init_distr_f, &param_f, &lum_model, &list_lum, &decay,
		   &param_decay, &time_of_run, &birthrate);
    }
    else                {
	in_param.open ("parameters.par");
	input_syntax (&in_param, &star_distr, &arms, &init_distr_p, &param_p, 
                   &init_distr_f, &param_f, &lum_model, &list_lum, &decay,
		   &param_decay, &time_of_run, &birthrate);
    }

in_param.close();

print_param(&cout, star_distr, arms, init_distr_p, param_p, 
                   init_distr_f, param_f, lum_model, list_lum, decay,
		   param_decay, time_of_run, birthrate);
	

    double T = -abs(time_of_run);    //-350e6;              // начало
    int number_stars = birthrate;    // темп звездообразования (звёзд в тысячалетие)
    int number_millenium = -T/1e3;   // количество тысячалетий


MFD * decay_model;
LM  * lum_mod;
RDF * rad_distr;
GD  * p_distr;
GD  * b_distr;

    if (type_of_run != 5)							{
	switch (star_distr)	{
		case 'A': rad_distr = new RDFFaucher; break;
 		case 'B': rad_distr = new RDFKruit;   break;
		case 'C': rad_distr = new RDFB0;      break;
		case 'D': rad_distr = new RDFSN;      break;
		case 'E': rad_distr = new RDFPuls;    break;
		default : print_error_no_RDF(star_distr); return 2; break;
	};

    rad_distr->print_description (&cout);

    }


    if (decay == 'A')	
    decay_model = new MFDConst (&param_decay);
    else if (decay == 'B')
    decay_model = new MFDStep  (&param_decay);
    else if (decay == 'C')
    decay_model = new MFDPons  (&param_decay);
    else if (decay == 'D')
    decay_model = new MFDExpon (&param_decay);
    else	{
    print_error_no_MFD(decay);
    return 6;
    }

    if (lum_model == 'B')
    lum_mod = new LMExpon (&list_lum);
    else if (lum_model == 'C')
    lum_mod = new LMFlat  (&list_lum);
    else	{
    print_error_no_LM(lum_model);	
    return 7;
    }
//    parametrs_lum param_lum ();
//    else



    decay_model->print_description  (&cout);
    decay_model->print_parameters    (&cout);
    lum_mod->print_description(&cout);
    lum_mod->print_parameters  (&cout);

    if (init_distr_p == 'g')
	p_distr = new GDGauss (&param_p);
    else if (init_distr_p == 'm')
	p_distr = new GDMGauss (&param_p);
    else  {
	print_error_flag_non_recognised('P');
	return 5;
   }

    if (init_distr_f == 'g')
	b_distr = new GDGauss (&param_f);
    else if (init_distr_f == 'm')
	b_distr = new GDMGauss (&param_f);
    else  {
	print_error_flag_non_recognised('B');
	return 5;
   }

    cout<<"#Distribution of initial periods."<<endl;
    p_distr->print_param        (&cout);
    cout<<"#Distribution of initial magnetic fields."<<endl;
    b_distr->print_param        (&cout);

    cout<<"#// Parameters of synthesis:                                    "<<endl;
    cout<<"#// T_start "<<T<<endl;
    cout<<"#// star formation rate "<<number_stars<<endl;
    cout<<"#//----------------------------------------------------------//"<<endl;
    cout<<"#// We are starting computations                             //"<<endl;
    srand(time(0));
    TMap T_copy;

    cout<<"#// Computations have been started                           //"<<endl;
    cout<<"#//----------------------------------------------------------//"<<endl;

ofstream out_res;
ofstream out_res_p;
ofstream out_res_pos;
ofstream out_res_lum;
ofstream out_res_hid;
ofstream out_res_p_hid;
ofstream out_res_pos_hid;
ofstream out_res_lum_hid;
ofstream out_bh;

char file_res_p[40], file_res_pos[40], file_res_lum[40];
char file_res_p_hid[40], file_res_pos_hid[40], file_res_lum_hid[40];
char file_res_hid[40];

    if (file_for_print)						{

	sprintf (file_res_p,   "p_%s"  , file_res);
	sprintf (file_res_pos, "pos_%s", file_res);
	sprintf (file_res_lum, "lum_%s", file_res);
	sprintf (file_res_hid,   "hidden_catalog_%s"  , file_res);
	sprintf (file_res_p_hid,   "hidden_catalog_p_%s"  , file_res);
	sprintf (file_res_pos_hid, "hidden_catalog_pos_%s", file_res);
	sprintf (file_res_lum_hid, "hidden_catalog_lum_%s", file_res);

	cout<<"#Result are written in "<<file_res<<endl;
	cout<<"#Result are written in "<<file_res_p<<endl;
	cout<<"#Result are written in "<<file_res_pos<<endl;
	cout<<"#Result are written in "<<file_res_lum<<endl;
   }

    if (!ext_print)
	if (file_for_print)
	    out_res.open (file_res);
        else 
	    out_res.open ("result.txt");
    else	{
	if (file_for_print)	{
	    out_res.open     (file_res);
	    out_res_p.open   (file_res_p);
	    out_res_pos.open (file_res_pos);
	    out_res_lum.open (file_res_lum);
	}
        else			{
	    out_res.open     ("result.txt");
	    out_res_p.open   ("p_result.txt");
	    out_res_pos.open ("pos_result.txt");
	    out_res_lum.open ("lum_result.txt");
	}
    }

    if (ext_print && hidden_file && file_for_print)	{
	    out_res_hid.open     (file_res_hid);
	    out_res_p_hid.open   (file_res_p_hid);
	    out_res_pos_hid.open (file_res_pos_hid);
	    out_res_lum_hid.open (file_res_lum_hid);	
    }
    else if (hidden_file && ext_print)			{
	    out_res_hid.open     ("hidden_catalog_result.txt");
	    out_res_p_hid.open   ("hidden_catalog_p_result.txt");
	    out_res_pos_hid.open ("hidden_catalog_pos_result.txt");
	    out_res_lum_hid.open ("hidden_catalog_lum_result.txt");
    }
    else if (hidden_file)
	    out_res_hid.open     ("hidden_catalog_result.txt");
	
    if (black_hole)
	    out_bh.open (BH_file);



    SpecialStar sun;
    SpecialStar sun_nowaday;
    double now = 0, shift;
    double P, dot_P, x, y, z, v_x, v_y, v_z, B, m, t1, t2;
    double dist_to_sun, lumin;
    double x_ances, y_ances, z_ances;
    int counter = 0;
    int rand_shift;                 // случайное время рождения внутри тысячелетия
    int n_magnet = 0;
    sun_nowaday.move_to(now);

    print_head (&out_res);

    out_res<<"#// Parameters of model:                                        "<<endl;
    out_res<<"#// T_start "<<T<<endl;
    out_res<<"#// star formation rate "<<number_stars<<endl;
    decay_model->print_description  (&out_res);
    decay_model->print_parameters    (&out_res);
    lum_mod->print_description(&out_res);
    lum_mod->print_parameters  (&out_res);
    p_distr->print_param        (&out_res);
    b_distr->print_param        (&out_res);
    out_res<<"#//----------------------------------------------------------//"<<endl;
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    out_res<<"#// Synthesis is started at "<<asctime(timeinfo);


    OBStar      * ancester;
    NeutronStar * descendant;

if (type_of_run == 4)					{

ofstream out_res (argum);

    for (int i = 0; i < number_millenium; i++)		{
        sun.move_to(T);

        for (int j = 0; j < number_stars; j++)	{
            rand_shift = rand()%500;
            ancester = new OBStar (T + rand_shift, &sun, rad_distr, arms);

            if (!(i%10000) && j==0) {
                cout<<T<<"\t"<<i<<endl;
            }

    shift  = ancester->get_time_on_MS();
    shift += ancester->get_t_He      ();

    ancester->move_to (shift);

    x = ancester->get_position_x ();
    y = ancester->get_position_y ();
    z = ancester->get_position_z ();

    v_x = ancester->get_velocity_x ();
    v_y = ancester->get_velocity_y ();
    v_z = ancester->get_velocity_z ();

    m = 1.17 + 0.09 * ancester->get_M_c_SN();

    t1 = T + rand_shift + shift;
    t2 = 0;

    out_res<<m<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<v_x/lcm*lsec/1e5<<"\t"<<v_y/lcm*lsec/1e5<<"\t"<<v_z/lcm*lsec/1e5;
    out_res<<t1<<"\t"<<t2<<endl;
   
    }

T += 1000;

}
return 0;
}
else if (type_of_run == 5)	{

   if (ext_print)	{	
	out_res_p   << "# P (s) \t dotP (s/s)"<<endl; 
	out_res_pos << "# x (kpc) \t y (kpc) \t z (kpc) \t v_x (km/s) \t v_y (km/s)\t v_z (km/s)"<<endl;
	out_res_lum << "# Lum at 1400 MHz (mJy kpc^2)"<<endl;
   }
   else
	out_res << "# P (s) \t dotP (s/s) \t x (kpc) \t y (kpc) \t z (kpc) \t v_x (km/s) \t v_y (km/s)\t v_z (km/s) \t Lum at 1400 MHz (mJy kpc^2)"<<endl;

cout<<"#Read file: "<<read_file<<endl;

ifstream in;
in.open (read_file);


int shet_i = 0;

    do { 
	in>>m;
	in>>x;
	in>>y;
	in>>z;
	in>>v_x;
	in>>v_y;
	in>>v_z;
	in>>t1;
	in>>t2;

	v_x = v_x / lsec*1e5*lcm;
 	v_y = v_y / lsec*1e5*lcm;
	v_z = v_z / lsec*1e5*lcm;

//cout<<shet_i<<endl;
	shet_i++;

	if (!(shet_i%10000))
		cout<<"line: "<<shet_i<<endl;

            descendant = new NeutronStar (t1-t2, decay_model, lum_mod, p_distr, b_distr, m, x, y, z, v_x, v_y, v_z, t2);
            P     = descendant->get_P(now);
//            out_str << descendant->get_velocity_x()<<endl;

		// if black_hole == true then write mass of BH in file
		if (!(descendant->is_this_ns()) && black_hole)
			out_bh << descendant->get_M()<<endl;

		// else: standard pulsar routine
                if (descendant->is_pulsar_alive(now) && descendant->is_this_ns() && lum_mod->is_beam_on(P))	{
			
                    lumin = descendant->is_pulsar_visible(now, &sun_nowaday, &T_copy);

                    if (lumin)				{
                        counter++;
                        x     = descendant->get_position_x();
                        y     = descendant->get_position_y();
                        z     = descendant->get_position_z();
                        v_x   = descendant->get_velocity_x();
                        v_y   = descendant->get_velocity_y();
                        v_z   = descendant->get_velocity_z();
                        P     = descendant->get_P(now);
                        dot_P = descendant->get_dot_P(now);
                        B     = descendant->get_B (now);
                        dist_to_sun = descendant->get_dist_to_sun(now, &sun_nowaday);

			if (ext_print)	{
        	                out_res_p  << P <<"\t"<< dot_P<<endl;
	                        out_res_pos<< x << "\t" << y <<"\t"<< z <<"\t"<<v_x<<"\t"<<v_y<<"\t"<<v_z<<endl;
	                        out_res_lum<<lumin*1000*pow(dist_to_sun, 2)<<endl;
				
				}
			else {
        	                out_res << P <<"\t"<< dot_P<<"\t";
	                        out_res << x << "\t" << y <<"\t"<< z <<"\t"<<v_x<<"\t"<<v_y<<"\t"<<v_z<<"\t";
	                        out_res <<lumin*1000*pow(dist_to_sun, 2)<<endl;
			}

                        if (B > 1e14) {
                            n_magnet++;
                        }
                    }
		   else if (hidden_file)		{
			x     = descendant->get_position_x();
			y     = descendant->get_position_y();
			z     = descendant->get_position_z();
			v_x   = descendant->get_velocity_x();					
			v_y   = descendant->get_velocity_y();
			v_z   = descendant->get_velocity_z();
			P     = descendant->get_P(now);
			dot_P = descendant->get_dot_P(now);
			B     = descendant->get_B (now);
			dist_to_sun = descendant->get_dist_to_sun(now, &sun_nowaday);

				if (ext_print)	{
			        out_res_p_hid  << P <<"\t"<< dot_P<<endl;
	                        out_res_pos_hid<< x << "\t" << y <<"\t"<< z <<"\t"<< v_x<<"\t"<<v_y<<"\t"<<v_z<<endl;
	                        out_res_lum_hid<<lumin*1000*pow(dist_to_sun, 2)<<endl;
				}
				else {
			        out_res_hid << P <<"\t"<< dot_P<<"\t";
	                        out_res_hid << x << "\t" << y <<"\t"<< z <<"\t"<< v_x<<"\t"<<v_y<<"\t"<<v_z<<"\t";
	                        out_res_hid <<lumin*1000*pow(dist_to_sun, 2)<<endl;
				}
			}

                }

                delete descendant;

     } while (!in.eof());

    time (&rawtime);
    timeinfo = localtime(&rawtime);
    out_res<<"#// Finished at "<<asctime(timeinfo);
    out_res<<"#// During synthesis it was detected as many pulsars   as  "<<counter<<endl;
    out_res<<"#// During synthesis it was detected as many magnetars as "<<n_magnet<<endl;

} 
else if (type_of_run == 2) {

   if (ext_print)	{	
	out_res_p   << "# P (s) \t dotP (s/s)"<<endl; 
	out_res_pos << "# x (kpc) \t y (kpc) \t z (kpc) \t v_x (km/s) \t v_y (km/s)\t v_z (km/s)"<<endl;
	out_res_lum << "# Lum at 1400 MHz (mJy kpc^2)"<<endl;
   }
   else
	out_res << "# P (s) \t dotP (s/s) \t x (kpc) \t y (kpc) \t z (kpc) \t v_x (km/s) \t v_y (km/s)\t v_z (km/s) \t Lum at 1400 MHz (mJy kpc^2)"<<endl;

	for (int i = 0; i < number_millenium; i++)		{
		sun.move_to(T);
		for (int j = 0; j < number_stars; j++)	{
			rand_shift = rand()%500;
			ancester = new OBStar (T + rand_shift, &sun, rad_distr, arms);

			if (!(i%2000) && j==0)
				cout<<T<<"\t"<<i<<endl;

			shift  = ancester->get_time_on_MS();
			shift += ancester->get_t_He      ();
			
			if (shift + T + rand_shift < 0)		{
				ancester->move_to (shift); // На сколько нужно сдвинуть, действительно мы не знаем времени рождения звезды
				descendant = new NeutronStar (T + shift + rand_shift, ancester, decay_model, lum_mod, p_distr, b_distr);
				P     = descendant->get_P(now);

				// if we need BH mass then print them
				if (!(descendant->is_this_ns()) && black_hole && i<10000)
				out_bh << descendant->get_M()<<endl;

				// else - standard pulsar routine
				if (descendant->is_pulsar_alive(now) && descendant->is_this_ns() && lum_mod->is_beam_on(P))	{
					descendant->move_to(now);
					lumin = descendant->is_pulsar_visible(now, &sun_nowaday, &T_copy);


					if (lumin)				{
						counter++;
						x     = descendant->get_position_x();
						y     = descendant->get_position_y();
						z     = descendant->get_position_z();
						v_x   = descendant->get_velocity_x();					
						v_y   = descendant->get_velocity_y();
						v_z   = descendant->get_velocity_z();
						P     = descendant->get_P(now);
						dot_P = descendant->get_dot_P(now);
						B     = descendant->get_B (now);
						dist_to_sun = descendant->get_dist_to_sun(now, &sun_nowaday);

						if (ext_print)	{
        	        			        out_res_p  << P <<"\t"<< dot_P<<endl;
				                        out_res_pos<< x << "\t" << y <<"\t"<< z <<"\t"<< v_x<<"\t"<<v_y<<"\t"<<v_z<<endl;
				                        out_res_lum<<lumin*1000*pow(dist_to_sun, 2)<<endl;
				
						}
						else {
        	        			        out_res << P <<"\t"<< dot_P<<"\t";
				                        out_res << x << "\t" << y <<"\t"<< z <<"\t"<< v_x<<"\t"<<v_y<<"\t"<<v_z<<"\t";
				                        out_res <<lumin*1000*pow(dist_to_sun, 2)<<endl;
						}

						if (B > 1e14)
							n_magnet++;		
						}
					   else if (hidden_file)		{
						x     = descendant->get_position_x();
						y     = descendant->get_position_y();
						z     = descendant->get_position_z();
						v_x   = descendant->get_velocity_x();					
						v_y   = descendant->get_velocity_y();
						v_z   = descendant->get_velocity_z();
						P     = descendant->get_P(now);
						dot_P = descendant->get_dot_P(now);
						B     = descendant->get_B (now);
						dist_to_sun = descendant->get_dist_to_sun(now, &sun_nowaday);

							if (ext_print)	{
        	        			        out_res_p_hid  << P <<"\t"<< dot_P<<endl;
				                        out_res_pos_hid<< x << "\t" << y <<"\t"<< z <<"\t"<< v_x<<"\t"<<v_y<<"\t"<<v_z<<endl;
				                        out_res_lum_hid<<lumin*1000*pow(dist_to_sun, 2)<<endl;
				
							}
							else {
        	        			        out_res_hid << P <<"\t"<< dot_P<<"\t";
				                        out_res_hid << x << "\t" << y <<"\t"<< z <<"\t"<< v_x<<"\t"<<v_y<<"\t"<<v_z<<"\t";
				                        out_res_hid <<lumin*1000*pow(dist_to_sun, 2)<<endl;
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
out_res<<"#// Finished at "<<asctime(timeinfo);
out_res<<"#// During synthesis it was detected as many pulsars   as  "<<counter<<endl;
out_res<<"#// During synthesis it was detected as many magnetars as "<<n_magnet<<endl;




}
    return 0;
}
