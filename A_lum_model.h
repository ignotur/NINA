#ifndef A_lum_model
#define A_lum_model

#include <fstream>
#include <iostream>

using namespace std;

class parametrs_lum {
public:
double ds, dlum;

parametrs_lum (ifstream *);

void print_description (ostream *);
void print_parametrs   (ostream *);
void print_short       (ostream *);

double get_ds   ();
double get_dlum (); 
 
private:

};

#endif
