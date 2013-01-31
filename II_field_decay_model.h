#ifndef II_field_decay_model
#define II_field_decay_model

#include <fstream>
#include <iostream>

using namespace std;

class parametrs_B {
public:
    double time_I, time_II, step;

    parametrs_B (ifstream *);

    void print_description (ostream *);
    void print_parametrs   (ostream *);
    void print_short       (ostream *);

    double get_time_I (void);
    double get_time_II(void);
    double get_step   (void);

private:

};

#endif
