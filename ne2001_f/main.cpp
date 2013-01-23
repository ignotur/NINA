#include <iostream>
#include <cmath>

using namespace std;

extern "C" {
void dmdsm_ (float *l, float *b, int *ndir, float *dmpsr, float *dist, char *limit, float *sm, float *smtau, float *smtheta, float *smiso);
}

int main () {
float pi = 3.1415926;
float l = 108.17/180.*pi, b = -42.98/180*pi, dist = 0.7; 
int ndir = -5;

float dmpsr, sm, smtau, smtheta, smiso;

char limit;
limit = ' ';

dmdsm_ (&l, &b, &ndir, &dmpsr, &dist, &limit, &sm, &smtau, &smtheta, &smiso);

cout<<dmpsr*sin(b)<<endl;

return 0;
}
