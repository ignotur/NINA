#ifndef IV_field_decay_model
#define IV_field_decay_model

#include <fstream>
#include <iostream>

using namespace std;

class parametrs_B {
public:
    double tau_ohm;

    parametrs_B (ifstream *);
    //parametrs_B (parametrs_B *);

    void print_description (ostream *);
    void print_parametrs   (ostream *);
    void print_short       (ostream *);

    //double get_alpha   ();
    double get_tau_ohm ();

private:

};

#endif
