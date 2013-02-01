#ifndef I_field_decay_model
#define I_field_decay_model

#include <fstream>
#include <iostream>

using namespace std;

class parametrs_B {
public:

    parametrs_B (ifstream *);

    void print_description (ostream *);
    void print_parametrs   (ostream *);
    void print_short       (ostream *);

private:

};

#endif
