#include <cmath>
#include "stars.h"

using namespace std;

double NeutronStar::get_P (double T)		{
double res;
res = base_mfd->get_P(T-tau, B, i_incl, P);
return res;
}

double NeutronStar::get_dot_P (double T)	{
double res;
res = base_mfd->get_dot_P(T-tau, B, i_incl, P);
return res;
}

double NeutronStar::get_B    (double T)		{
double res;
res = base_mfd->get_B(T-tau, B, i_incl, P);
return res;
}

double NeutronStar::is_pulsar_visible (double t, SpecialStar * sun, TMap * T_copy) {
double res, P_curr, dot_P_curr, DM, dist_to_sun;
double first [3], second[2];
float l, b, sm;
P_curr     = base_mfd->get_P    (t-tau, B, i_incl, P);
dot_P_curr = base_mfd->get_dot_P(t-tau, B, i_incl, P);

    sun->move_to(t);
    dist_to_sun = sqrt(pow(sun->get_position_x() - x, 2) + pow(sun->get_position_y() - y, 2) + pow(sun->get_position_z() - z, 2));

    first [0] = x - sun->get_position_x();
    first [1] = y - sun->get_position_y();
    first [2] = z - sun->get_position_z();


    // Вектор от Солнца к центру Галактики
    second[0] = - sun->get_position_x();
    second[1] = - sun->get_position_y();

    b = first[2] / dist_to_sun;
    b = asin (b)/2./pi*360.;

    l = atan2(second[0]*first[1]-second[1]*first[0], first[0]*second[0]+first[1]*second[1])*180./3.1415926;

    if (l<0) {
        l=360.+l;
    }

DM         = get_DM (t, sun, &l, &b, &sm);

res = base_lm->is_pulsar_visible (t-tau, sun, T_copy, x, y, z, i_incl, P_curr, dot_P_curr, DM);
return res;
}
