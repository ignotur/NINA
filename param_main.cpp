#include <iostream>
#include <fstream>
#include "stars.h"


using namespace std;

int main () {

ifstream input ("input.txt");

parametrs_B    param_B   (&input);
parametrs_lum  param_lum (&input);

param_B.print_description  (&cout);
param_B.print_parametrs    (&cout);
param_lum.print_description(&cout);
param_lum.print_parametrs  (&cout);

return 0;
}
