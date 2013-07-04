#include <iostream>
#include <vector>
using namespace std;

double version () { 
return 0.8;
}

void print_head (ostream * out)	{
*out<<"#//----------------------------------------------------------//"<<endl;
*out<<"#// Пульсарная популяция. Версия "<<version()<<" см соглашение о версиях //"<<endl;
*out<<"#// Автор: Игошев Андрей, научный руководитель: А.Ф. Холтыгин//"<<endl;
*out<<"#// e-mail: ignotur@gmail.com СПбГУ, 2010-2013               //"<<endl;
*out<<"#//----------------------------------------------------------//"<<endl;
}

void print_help ()	{
cout<<"The population synthesis code allows to use following flags    "<<endl;
cout<<"-e -- extended output: one output file is divided into 3 files."<<endl;
cout<<"      One with P-dotP, one with positions and velocities and   "<<endl;
cout<<"      one with luminosities at 1400 MHz.                       "<<endl;
cout<<"-f -- full population. These type of calculations starts from  "<<endl;
cout<<"      massive stars in four spiral arms, then stars evolve and "<<endl;
cout<<"      produce NS. Their evolution is considered up to now      "<<endl;
cout<<"-g file -- generate file with population of massive stars.     "<<endl;
cout<<"           This file constain masses of newly-born NS, their   "<<endl;
cout<<"           positions, velocities and time of born              "<<endl;
cout<<"-h -- help - show this information.                            "<<endl;
cout<<"-i file -- initial file. This flag makes program use file with "<<endl;
cout<<"           initial positions, velocities and masses of NS.     "<<endl;
cout<<"-o file -- file/s with results has/ve name as suggested.       "<<endl;
cout<<"-p file -- file with parameters of synthesis is names as       "<<endl;
cout<<"           suggested. If it is not mentioned, file input.par is"<<endl;
cout<<"           used.                                               "<<endl;
cout<<"-u -- undetectable. This option makes program create additional"<<endl;
cout<<"      file/s with information about pulsars missed in survey.  "<<endl;
cout<<"      (This option may cause much longer computations)         "<<endl;
}

void print_exclusion () {
cout<<"Error 1: flags f, g and i cannot be used simultaneously."<<endl;
}

void print_param (ostream * out, char star_distr, bool arms, char init_distr_p, vector<double> param_p, 
      char init_distr_f, vector<double> param_f, char lum_model, vector <double> param_lum, char decay,
		   vector<double> param_decay, double time_of_run, double birthrate_val)  {

*out<<"---------------------------------------------------------------"<<endl;
*out<<"Parameters of run."<<endl;
*out<<"Interval of simulation "<<time_of_run<<" years."<<endl;
*out<<"Birthrate of massive stars "<<birthrate_val<<" per thousand of years."<<endl;
*out<<"Model of initial radial distribution of stars "<<star_distr<<endl;
*out<<"Initial distribution of pulsars periods "<<init_distr_p<<endl;
*out<<"Initial distribution of pulsars's magnetic fields "<<init_distr_f<<endl;
*out<<"Model of luminocity "<<lum_model<<endl;
*out<<"Model of magnetic field decay "<<decay<<endl;
}

void print_error_no_RDF(char tmp) {
cout<<"Error 2: no model of with identifier "<<tmp<<" was specified. Check file with parameters."<<endl;
}

void print_error_parameters_not_enough () {
cout<<"Error 3: not enough parameters are specified for generating of initial distributions."<<endl;
}

void print_error_flag_non_recognised (char val) {
cout<<"Error 5: in file with parameters flag for the distribution of initial "<<val<<" is not recognised."<<endl;
}
