#ifndef C_lum_model
#define C_lum_model

#include <fstream>
#include <iostream>

using namespace std;

class parametrs_lum {
public:

parametrs_lum (ifstream *);

bool is_beam_on(double);

void print_description (ostream *);
void print_parametrs   (ostream *);
void print_short       (ostream *);
 
private:

};

#endif
