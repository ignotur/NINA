#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <functional>
#include "stars.h"

using namespace std;

//static_cast<int(*)(int)>(isspace);

void input_syntax (ifstream * in,  char * star_distr, bool * arms, char * init_distr_p, vector<double>* param_p, 
      char * init_distr_f, vector<double>* param_f, char * lum_model, vector <double>* param_lum, char * decay,
		   vector<double>* param_decay, double * time_of_run, double * birthrate_val)	{

string data[100];
string A1;
int n=0;

	do {
		getline(*in, data[n], ':');
//cout<<data[n]<<endl;
		n++;
		getline(*in, data[n]);
//cout<<data[n]<<endl;
		n++;
	} while (!in->eof());
n--; n--;

//cout<<"----------------------------------------------------------------"<<endl;

string initial   = "Initial distribution";
string MFd       = "Model of magnetic field decay";
string radial_star_distr = "Initial radial distribution";
string lum       = "Model of luminocity";
string birthrate = "Birthrate";
string time_run  = "Time of run";
string param_per = "parameter of initial distribution of pulsars periods";
string period    = "Initial distribution of pulsars period";
string spiral    = "Spiral arms";
string param_mag = "parameter of initial distribution of pulsars magnetic fields";
string field     = "Initial distribution of pulsars magnetic fields";
string param_l   = "parameter of luminosity";
string param_mfd = "parameter of magnetic field decay";

string str;

const char * res_compar;
const char * tmp, *tmp_comp;
const char yes[] = "yes";

	for (int i=0; i<n/2; i++)	{
	tmp      = data[i*2].c_str();

	tmp_comp = time_run.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
		*time_of_run = atof(data[i*2+1].c_str());
//	cout<<data[i*2]<<"Time of run string"<<endl;	
//	cout<<*time_of_run<<endl;
	
	}

	tmp_comp = birthrate.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
		*birthrate_val = atof(data[i*2+1].c_str());
//	cout<<data[i*2]<<"Time of run string"<<endl;	
//	cout<<*birthrate_val<<endl;
	
	}

	tmp_comp = radial_star_distr.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
	 	data[i*2+1].erase(remove_if(data[i*2+1].begin(), data[i*2+1].end(), ptr_fun <int, int> ( isspace ) ), data[2*i+1].end());
		*star_distr = data[i*2+1].at(0);
//	cout<<data[i*2]<<"Time of run string"<<endl;	
//	cout<<*star_distr<<endl;
	
	}

	tmp_comp = lum.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
	 	data[i*2+1].erase(remove_if(data[i*2+1].begin(), data[i*2+1].end(), ptr_fun <int, int> ( isspace ) ), data[2*i+1].end());
		*lum_model = data[i*2+1].at(0);
//	cout<<data[i*2]<<"Time of run string"<<endl;	
//	cout<<*star_distr<<endl;
	
	}

	tmp_comp = MFd.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
	 	data[i*2+1].erase(remove_if(data[i*2+1].begin(), data[i*2+1].end(), ptr_fun <int, int> ( isspace ) ), data[2*i+1].end());
		*decay = data[i*2+1].at(0);
//	cout<<data[i*2]<<"Time of run string"<<endl;	
//	cout<<*star_distr<<endl;
	
	}

	tmp_comp = spiral.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
		tmp_comp = data[2*i+1].c_str();
		res_compar = strstr(tmp_comp, yes);
		if (res_compar != NULL)
			*arms = true;
		else
			*arms = false;
//		cout<<*arms<<endl;	
	}

	tmp_comp = period.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
	 	data[i*2+1].erase(remove_if(data[i*2+1].begin(), data[i*2+1].end(), ptr_fun <int, int> ( isspace ) ), data[2*i+1].end());
		*init_distr_p = data[i*2+1].at(0);
	}

	tmp_comp = param_per.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
		param_p->push_back(atof(&data[2*i].at(0)));
		param_p->push_back(atof(data[i*2+1].c_str()));

//cout<<"Here"<<endl;
//cout<<atof(&data[2*i].at(0)) << "\t" << atof(data[i*2+1].c_str())<<endl;

	}

	tmp_comp = field.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
	 	data[i*2+1].erase(remove_if(data[i*2+1].begin(), data[i*2+1].end(), ptr_fun <int, int> ( isspace ) ), data[2*i+1].end());
		*init_distr_f = data[i*2+1].at(0);
	}

	tmp_comp = param_mag.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
		param_f->push_back(atof(&data[2*i].at(0)));
		param_f->push_back(atof(data[i*2+1].c_str()));
	}

	tmp_comp = param_l.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
		param_lum->push_back(atof(&data[2*i].at(0)));
		param_lum->push_back(atof(data[i*2+1].c_str()));
	}

	tmp_comp = param_mfd.c_str(); 
	res_compar = strstr(tmp, tmp_comp);
	if (res_compar != NULL)					{
		param_decay->push_back(atof(&data[2*i].at(0)));
		param_decay->push_back(atof(data[i*2+1].c_str()));
	}

//	tmp_comp = initial.c_str(); 
//	res_compar = strstr(tmp, tmp_comp);
//	if (res_compar != NULL)					{
//	cout<<data[i*2]<<" is initial string"<<endl;	
//	}

	}

}
