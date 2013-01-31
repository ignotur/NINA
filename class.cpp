#include <cstdlib>
#include <time.h>
#include <cmath>
#include <iostream>
#include "stars.h"
#include <algorithm>

using namespace std;

//------------------------------------------------------------------------------//
// star_OB --- класс, который позволяет создать и перемешать массивную звезды в
// Галактике. Конструктор класса создает звезду расположенную в окрестности
// спирального рукава, с параметрами генерируемыми случайно, но в заданных
// физически обоснованных  пределах
//
//  Автор: Igoshev Andrey
//  Научный руководитель: Холтыгин А.Ф.
//  e-mail: igoshev-andrei@rambler.ru
//  Написание начато: 13.09.2010
//------------------------------------------------------------------------------//


//------------------------------------------------------------------------------//
// Замечание об используемых единицах измерения:
// Единица массы      - грамм
// Единица времени    - год = 365.24*24*3600 сек
// Единица расстояния - кпк = 3.085678e+21 см
//------------------------------------------------------------------------------//


//--------------------------------------//
// Декларация функций
double rho (double);
double rho_P90 (double);
double rho_B00 (double);
double rho_SN_remnant(double);
double rho_F06 (double);
double rho_m (double);
double norm_distr (void);
double dphi_dx (double, double, double);
double dphi_dy (double, double, double);
double dphi_dz (double, double, double);
double dphi_dr (double, double, double);
void   Runge_Kutta (int, double, double *, void (*f)(int, double *));
void   diff_equi   (int, double *);
void   set_position(double, double, double);
void   set_velocity(double, double, double);
//--------------------------------------//

OBStar::OBStar (double T, SpecialStar * sun) {
    double chance_1, chance_2, theta_corr, theta, corr; //theta - угол из цетра Галактики
    int    arm;
    bool is_position_set = false;
    double r, dx, dy, s;
    double v;

    // Параметры спиральных рукавов Галактики ---
    // в настоящий момент времени             ---
    double k[4], r_0[4], theta_0 [4];


    //-----------------------------------------------//
    //k[0] = 4.25; r_0[0] = 3.48; theta_0[0] = 1.57; // Norma
    //k[1] = 4.25; r_0[1] = 3.48; theta_0[1] = 4.71; // Carina-Sagittarius
    //k[2] = 4.89; r_0[2] = 4.90; theta_0[2] = 4.09; // Perseus
    //k[3] = 4.89; r_0[3] = 4.90; theta_0[3] = 0.95; // Crux-Scutum
    // --------------------------------------------- //

    k[0] = 5.33;
    r_0[0] = 2.52;
    theta_0[0] = pi/2.; // Norma
    k[1] = 5.33;
    r_0[1] = 2.52;
    theta_0[1] = pi; // Carina-Sagittarius
    k[2] = 5.33;
    r_0[2] = 2.52;
    theta_0[2] = 3*pi/2.; // Perseus
    k[3] = 5.33;
    r_0[3] = 2.52;
    theta_0[3] = 4*pi/2.; // Crux-Scutum

    arm = rand()%4;

    // Генерация радиального распределения в соотвествии со статьёй
    // Faucher-Giguere, Kaspi, 2006 методом усечения.
    // Диаметр Галактики взят 30 кпк

    do {

        x = rand () / rand_high_board;
        x*=12.48;
        x+=2.52;
        y = rand () / rand_high_board;

        if (rho(x) >= y) {            ///rho_P90
            is_position_set = true;
        }

    } while (!(is_position_set));

    //x-=2.5;
    //x*=3;
    r = x;

    theta = k[arm] * log(r/r_0[arm]) + theta_0[arm];
    theta_corr =  rand() / rand_high_board;
    theta_corr *= 2 * pi;
    sun->move_to(T);
    theta += theta_corr * pow(e, -3.5*r) + sun->get_theta();
    //theta=theta_corr;
    //---------------------------------------------------------------

    // Поворот спирального узора из-за вращения Галактики
    // используется статья Ramachandran, 1994
    // используется предположение твердотельного вращения спирального узора
    // с постоянной угловой скоростью равной 23 км/сек кпк
    // из работы Л.С. Марочник А.А. Сучков

    corr = 2.55673e-7/r*T;

    theta_corr = 2.301055e-7* (1./r - 1./14)*T; // это 225 км/с переведённые в кпк/год

    //********************************************************This is new
    theta_corr = 2.35219e-8*T;

    theta += theta_corr;

    //---------------------------------------------------------------

    x = r * cos(theta);
    y = r * sin(theta);

    // Генерация нормально распределённой случайной величины
    // методом Бокса-Мюллера для смазывания положения рождающихся
    // протегионов относительно центроидов спиральных рукавов

    is_position_set = false;

    do {
        chance_1 = rand () / rand_high_board;
        chance_2 = rand () / rand_high_board;
        chance_1 = 2*(chance_1 - 0.5);
        chance_2 = 2*(chance_2 - 0.5);
        s = chance_1*chance_1 + chance_2*chance_2;

        if (s != 0 && s<=1) {
            is_position_set = true;
        }
    } while (!(is_position_set));

    dx = chance_1 * sqrt(-2*log(s)/s);
    dy = chance_2 * sqrt(-2*log(s)/s);
    dx *= 0.07*r;
    dy *= 0.07*r;

    x+=dx;
    y+=dy;

    //---------------------------------------------------------------

    z = 0.050 * norm_distr ();

    // Генерация пространственных скоростей протегионов пульсаров
    // дисперсия во всех направлениях 15 км/сек из работы
    // Ramachandran, 1994; кривая вращения Галактики согласована
    // с её гравитационным потенциалом
    r = sqrt(x*x + y*y);
    v_z = 1.534e-8*norm_distr ();
    v = sqrt(r*dphi_dr(x,y,0));
    v_x = - v * sin(theta);
    v_y =   v * cos(theta);

    v_x += 1.534e-8*norm_distr();
    v_y += 1.534e-8*norm_distr();

    //----------------------------------------------------------------

    // Устанавливаем текущее время временем рождения + некоторый разброс внутри тысячялетия
    tau = 1000*rand() / rand_high_board;
    tau+= T;
    //----------------------------------------------------------------

    // Распределение магнитных потоков звёзд логарифм-нормальное,
    // по статье Игошев А.П. Холтыгин А.Ф., 2010

    F = 27.13 + 0.67 * norm_distr();
    F = pow(10,F);

    //----------------------------------------------------------------

    // Генерация звёзд с массой заключённой в пределах 8-45 M_sol
    // в соотвествии с её функцией плотности

    is_position_set = false;

    do {
        chance_1 = rand () / rand_high_board;
        chance_1 *= 37;
        chance_1 += 8;
        chance_2 = rand () / rand_high_board;

        if (rho_m(chance_1) >= chance_2) {
            is_position_set = true;
        }

    } while (!(is_position_set));

    M = chance_1;
    //----------------------------------------------------------------

    //----------------------------------------------------------------
    // Распределение рождающихся звёзд по металичностям
    // В окрестности солнца - 0.02, далее радиальный
    // галактический градиент как в статье Maciel & Costa 2010 года
    // по цефеидам, а именно: до 8 кпк -0.130+-0.015 кпк^-1
    // после 8 кпк -0.042+-0.004 кпк^-1

    r = sqrt(x*x+y*y);

    if (R<=8.) {
        Z = log10(0.02) - (8. - r)*(0.130 + 0.015 * norm_distr());
    } else {
        Z = log10(0.02) - (r - 8.)*(0.042 + 0.004 * norm_distr());
    }

    //Z = log10(0.002);
    Z = pow(10., Z);

    ksi = log10(Z/0.02);
    //----------------------------------------------------------------

    // Расчёт времени жизни звезды на главной последовательности и
    // времени пересечения Гершпрунгового пробела???

    //----------------------------------------------------------------

    // Расчёт радиуса звезды ZAMS исходя из аппроксимации выполненной в
    // статье Hurley, 2000

    double a17, sigma, c1, a[33][5] ;
    sigma = log10(Z);
    c1 = -8.672073e-2;
    a17 = max(0.097-0.1072*(sigma+3), max(0.097, min(0.1461, 0.1461+0.1237*(sigma+2))));

    a[23] [0] = 2.617890;
    a[23] [1] = 1.019135;
    a[23] [2] = -3.292551e-2;
    a[23] [3] = -7.445123e-2;
    a[23] [4] = 0;
    a[24] [0] = 1.075567e-2;
    a[24] [1] = 1.773287e-2;
    a[24] [2] = 9.610479e-3;
    a[24] [3] = 1.732469e-3;
    a[24] [4] = 0;
    a[25] [0] = 1.476246;
    a[25] [1] = 1.899331;
    a[25] [2] = 1.195010;
    a[25] [3] = 3.035051e-1;
    a[25] [4] = 0;
    a[26] [0] = 5.502535;
    a[26] [1] = -6.601663e-2;
    a[26] [2] = 9.968707e-2;
    a[26] [3] = 3.599801e-2;
    a[26] [4] = 0;

    for (int i = 23; i < 27; i++)
        for (int j = 1; j < 5; j++) {
            a[i][0] += a[i][j] * pow(ksi, j);
        }

    R = (c1*pow(M, 3) + a[23][0]*pow(M, a[26][0]) + a[24][0]*pow(M, a[26][0] + 1.5))/(a[25][0] + pow(M, 5));
    //R = log10(R);
    //--------------------------------------------------------------------
    // Расчёт угловой момента импульса
    // Скорость вращения звезды взята равной 200 км/с из К.У. Аллена  1987
    // В дальнейшем будем предполагать эту величину сохраняющейся со временем
    // Так же по статье Hurley, 2000
    //L = 2./5. * M*M_sol*R_sol*R*200e5;
    double v_rot, Omega;
    v_rot = 330*pow(M, 3.3)/(15.0+pow(M, 3.45));
    Omega = v_rot /R/R_sol*1e5;
    //Omega *= lsec;
    L = 2./5.*M*M_sol*R*R_sol*R*R_sol*Omega;
    //Omega = L / 2.*5./M/M_sol/pow(1000*R_sol, 2);

    //L = 4*pi*M*M_sol*pow(R*R_sol, 2)/5./L;
    //L = 2*pi/Omega;
    //---------------------------------------------------------------------

    // Расчёт времени жизни на Главной последовательности и времени прохождения
    // пробела Гершпрунга. По статье Hurley, 2000

    double hi, mu;
    double L_bgb_M_HeF, L_HeI_M_HeF;

    a[1] [0] = 1.593890e3;
    a[1] [1] = 2.053038e3;
    a[1] [2] = 1.231226e3;
    a[1] [3] = 2.327785e2;
    a[1] [4] = 0;
    a[2] [0] = 2.706708e3;
    a[2] [1] = 1.483131e3;
    a[2] [2] = 5.772723e2;
    a[2] [3] = 7.411230e1;
    a[2] [4] = 0;
    a[3] [0] = 1.466143e2;
    a[3] [1] = -1.04844e2;
    a[3] [2] = -6.795374e1;
    a[3] [3] = -1.391127e1;
    a[3] [4] = 0;
    a[4] [0] = 4.141960e-2;
    a[4] [1] = 4.564888e-2;
    a[4] [2] = 2.958542e-2;
    a[4] [3] = 5.571483e-3;
    a[4] [4] = 0;
    a[5] [0] = 3.426349e-1;
    a[5] [1] = 0;
    a[5] [2] = 0;
    a[5] [3] = 0;
    a[5] [4] = 0;
    a[6] [0] = 1.949814e1;
    a[6] [1] = 1.758178;
    a[6] [2] = -6.008212;
    a[6] [3] = -4.470533;
    a[6] [4] = 0;
    a[7] [0] = 4.903830;
    a[7] [1] = 0;
    a[7] [2] = 0;
    a[7] [3] = 0;
    a[7] [4] = 0;
    a[8] [0] = 5.212154e-2;
    a[8] [1] = 3.166411e-2;
    a[8] [2] = -2.750074e-3;
    a[8] [3] = -2.271549e-3;
    a[8] [4] = 0;
    a[9] [0] = 1.312179;
    a[9] [1] = -3.294936e-1;
    a[9] [2] = 9.231860e-2;
    a[9] [3] = 2.610989e-2;
    a[9] [4] = 0;
    a[10][0] = 8.073972e-1;
    a[10][1] = 0;
    a[10][2] = 0;
    a[10][3] = 0;
    a[10][4] = 0;

    for (int i = 1; i < 10; i++)
        for (int j = 1; j < 5; j++) {
            a[i][0] += a[i][j] * pow (ksi, j);
        }

    t_bgb = (a[1][0] + a[2][0]*pow(M, 4) + a[3][0] * pow(M, 5.5) + pow(M, 7))/(a[4][0] * pow(M, 2) + a[5][0] * pow(M, 7));
    hi = max(0.95, min(0.95 - 0.03*(ksi + 0.30103), 0.99));
    mu = max(0.5, 1.-0.01*max(a[6][0]/pow(M, a[7][0]), a[8][0] + a[9][0]/pow(M, a[10][0])));
    t_ms = max(mu * t_bgb, hi * t_bgb);
    //t_ms*=1e6;
    //t_bgb*=1e6;
    //---------------------------------------------------------------------------

    // Расчёт времени прохождения ветви гигантов, по статье Hurley, 2000
    double D, A_H, p = 5, q = 2, L_bgb, L_HeI, t_inf_1, t_x, t_inf_2, B, b17;
    double M_HeF, M_FGB, L_x, L_min_He, c;
    double c2 = 9.301992, c3 = 4.637345;
    D = max (max(-1.0, 0.975*(5.37+0.135*ksi)-0.18*M), 0.5*(5.37+0.135*ksi)-0.06*M);
    D = pow(10, D);
    A_H = max(-4.8, min(-5.7+0.8*M, -4.1+0.14*M));

    a[27][0] = 9.511033e+1;
    a[27][1] = 6.819618e+1;
    a[27][2] = -1.045625e+1;
    a[27][3] = -1.474939e+1;
    a[27][4] = 0;
    a[28][0] = 3.113458e+1;
    a[28][1] = 1.012033e+1;
    a[28][2] = -4.650511;
    a[28][3] = -2.463185;
    a[28][4] = 0;
    a[29][0] = 1.413057;
    a[29][1] = 4.578814e-1;
    a[29][2] = -6.850581e-2;
    a[29][3] = -5.588658e-2;
    a[29][4] = 0;
    a[30][0] = 3.910862e+1;
    a[30][1] = 5.196646e+1;
    a[30][2] =  2.264970e+1;
    a[30][3] =  2.873680;
    a[30][4] = 0;
    a[31][0] = 4.597479;
    a[31][1] =-2.855179e-1;
    a[31][2] =  2.709724e-1;
    a[31][3] =  0;
    a[31][4] = 0;
    a[32][0] = 6.682518;
    a[32][1] = 2.827718e-1;
    a[32][2] = -7.294429e-2;
    a[32][3] =  0;
    a[32][4] = 0;

    for (int i = 27; i < 33; i++)
        for (int j = 1; j < 5; j++) {
            a[i][0] += a[i][j] * pow (ksi, j);
        }

    a[29][0] = pow(a[29][0], a[32][0]);


    L_bgb = (a[27][0]*pow(M, a[31][0]) + a[28][0]*pow(M, c2))/(a[29][0] + a[30][0]*pow(M,c3) + pow(M, a[32][0]));

    t_inf = t_bgb + 1./(A_H*D*(p-1.)) * pow(D/L_bgb, (p-1.)/p);

    double b [45][4];

    b[11][0] = 1.071738e+2;
    b[11][1] = -8.970339e+1;
    b[11][2] = -3.949739e+1;
    b[12][0] = 7.348793e+2;
    b[12][1] = -1.531020e+2;
    b[12][2] = -3.793700e+1;
    b[13][0] = 9.219293e+0;
    b[13][1] = -2.005865e+0;
    b[13][2] = -5.561309e-1;
    b[14][0] = 2.917412e+0;
    b[14][1] =  1.575290e+0;
    b[14][2] =  5.751814e-1;
    b[15][0] = 3.629118e+0;
    b[15][1] = -9.112722e-1;
    b[15][2] =  1.042291e+0;
    b[16][0] = 4.916389e+0;
    b[16][1] =  2.862149e+0;
    b[16][2] =  7.844850e-1;

    for (int i = 11; i < 17; i++)
        for (int j = 1; j < 3; j++) {
            b[i][0] += b[i][j] * pow (ksi, j);
        }

    b[11][0] = pow(b[11][0], 2);
    b[13][0] = pow(b[13][0], 2);
    b[14][0] = pow(b[14][0], b[15][0]);
    b[16][0] = pow(b[16][0], b[15][0]);

    if (ksi>-1.0) {
        b17 = 1.-0.3880523*pow((ksi+1.0), 2.862149);
    } else {
        b17 = 1.;
    }

    M_HeF = 1.995+0.25*ksi + 0.087 * pow(ksi, 2);
    M_FGB = (13.048*pow(Z/0.02, 0.06))/(1+0.0012*pow(0.02/Z, 1.27));

    c = b17 / pow(M_FGB, 0.1) + (b[16][0]*b17 - b[14][0])/(pow(M_FGB, b[15][0]+0.1));

    L_HeI    = (b[11][0] + b[12][0]*pow(M, 3.8))/(b[13][0] + pow(M, 2));
    L_min_He = L_HeI * (b[14][0] + c * pow(M, b[15][0]+0.1))/(b[16][0] + pow(M, b[15][0]));


    B = max (3e+4, 500+1.75e+4*pow(M, 0.6));

    if (M > M_HeF && M < M_FGB) {
        L_x = L_min_He;
    } else if (M>=M_FGB) {
        L_x = L_HeI;
    } else {
        cout<<"Error, mass of star isn't in possible border. E.g. M = "<<M<<", M_HeF = "<<M_HeF<<", M_FGB = "<<endl;
    }

    t_inf_1 = t_bgb + 1./(p-1.) / A_H / D * pow(D/L_bgb, (p-1.)/p);
    t_x     = t_inf_1 - (t_inf_1 - t_bgb) * pow(L_bgb/L_x, (p-1.)/p);
    t_inf_2 = t_x + 1/(q-1) / A_H / B * pow(B/L_x, (q-1.)/q);

    if (L_x >= L_HeI) {
        t_HeI = t_inf_1 - 1./(p-1.)/A_H/D*pow(D/L_HeI, (p-1.)/p);
    } else if (L_HeI > L_x) {
        t_HeI = t_inf_2 - 1./(q-1.)/A_H/B*pow(B/L_HeI, (q-1.)/q);
    } else {
        cout<<"t_HeI wasn't set!"<<endl;
    }

    b[41][0] = 2.327037e+0;
    b[41][1] = 2.403445e+0;
    b[41][2] = 1.208407e+0;
    b[41][3] = 2.087263e-1;
    b[42][0] = 1.997378e+0;
    b[42][1] =-8.126205e-1;
    b[42][2] = 0;
    b[42][3] = 0;
    b[43][0] = 1.079113e-1;
    b[43][1] = 1.762409e-2;
    b[43][2] = 1.096601e-2;
    b[43][3] = 3.058818e-3;
    b[44][0] = 2.327409e+0;
    b[44][1] = 6.901582e-1;
    b[44][2] =-2.158431e-1;
    b[44][3] =-1.084117e-1;

    for (int i = 41; i < 45; i++)
        for (int j = 1; j < 4; j++) {
            b[i][0] += b[i][j]*pow(ksi, j);
        }


    b[41][0] = pow(b[41][0], b[42][0]);
    b[44][0] = pow(b[44][0], 5);

    t_He = t_bgb * (b[41][0]*pow(M, b[42][0]) + b[43][0]*pow(M,5))/(b[44][0] + pow(M, 5));

    // Расчёт массы ядра на разных стадиях эволюции массивной звезды

    L_bgb_M_HeF = (a[27][0]*pow(M_HeF, a[31][0]) + a[28][0]*pow(M_HeF, c2))/(a[29][0] + a[30][0]*pow(M_HeF,c3) + pow(M_HeF, a[32][0]));

    double M_Ch= 1.44, C;
    double GB3 = 4.1e+4;
    double GB4 = 5.5e+4/(1+0.4*pow(M,4));
    double GB5 = 5.;
    double GB6 = 3.;
    c1 = 9.20925e-5;
    c2 = 5.0402216;

    b[36][0] = 1.445216e-1;
    b[36][1] = -6.180219e-2;
    b[36][2] = 3.093878e-2;
    b[36][3] =  1.567090e-2;
    b[37][0] = 1.304129e+0;
    b[37][1] =  1.395919e-1;
    b[37][2] = 4.142455e-3;
    b[37][3] = -9.732503e-3;
    b[38][0] = 5.114149e-1;
    b[38][1] = -1.160850e-2;
    b[38][2] = 0;
    b[38][3] = 0;

    for (int i = 36; i < 39; i++)
        for (int j = 1; j < 4; j++) {
            b[i][0] += b[i][j]*pow(ksi, j);
        }

    b[36][0] = pow(b[36][0],4);
    b[37][0] = 4.*b[37][0];
    b[38][0] = pow(b[38][0],4);

    if (L_bgb_M_HeF < L_x) {
        M_c_HeF = pow(L_bgb_M_HeF/GB4, 1./GB5);
    } else {
        M_c_HeF = pow(L_bgb_M_HeF/GB3, 1./GB6);
    }

    M_c_BAGB = pow(b[36][0]*pow(M, b[37][0]) + b[38][0], 0.25);
    C = pow(M_c_HeF, 4) - c1*pow(M_HeF, c2);
    M_c_BGB  = min(0.95*M_c_BAGB, pow(C + c1*pow(M, c2), 0.25));

    L_HeI_M_HeF   = (b[11][0] + b[12][0]*pow(M_HeF, 3.8))/(b[13][0] + pow(M_HeF, 2));

    if (L_HeI_M_HeF < L_x) {
        M_c_HeI = pow(L_HeI_M_HeF/GB4, 1./GB5);
    } else {
        M_c_HeI = pow(L_HeI_M_HeF/GB3, 1./GB6);
    }

    C = pow(M_c_HeI, 4) - c1*pow(M_HeF, c2);
    M_c_HeI  = min(0.95*M_c_BAGB, pow(C + c1*pow(M, c2), 0.25));
    M_c_DU   = 0.44*M_c_BAGB + 0.448;
    M_c_SN   = max(M_Ch, 0.773*M_c_BAGB-0.35);

    //L -= Omega * (M-M_c_SN)*M_sol * 0.1 * pow(1000*R_sol, 2);


    //---------------------------------------------------------------------------
    // Потеря углового момента за счёт магнитного торможения, на поздних стадиях
    // эволюции звезды, когда образуеться атмосфера
    // Так же по статье Hurley, 2000
    /*double J = 0, M_env;
    M_env = M - M_c_BGB;
    J = 5.83e-16 * M_env/M*pow(R*Omega, 3);
    J *= t_bgb*1e6;
    L -= J;
    Omega = */
    //---------------------------------------------------------------------------
}



double OBStar::get_position_x () {
    return x;
}

double OBStar::get_position_y () {
    return y;
}

double OBStar::get_position_z () {
    return z;
}

double OBStar::get_velocity_x () {
    return v_x;
}

double OBStar::get_velocity_y () {
    return v_y;
}

double OBStar::get_velocity_z () {
    return v_z;
}

void OBStar::move_to(double T) {
    double result [6];

    result [0] = x;
    result [1] = y;
    result [2] = z;
    result [3] = v_x;
    result [4] = v_y;
    result [5] = v_z;

    //cout<<"A star OB is moving to "<<T<<endl;

    Runge_Kutta (6, T, &result[0], &diff_equi);

    x   = result [0];
    y   = result [1];
    z   = result [2];
    v_x = result [3];
    v_y = result [4];
    v_z = result [5];

}

double OBStar::get_Flux () {
    return F;
}

double OBStar::get_mass () {
    return M;
}

void OBStar::set_position (double x_, double y_, double z_) {
    x = x_;
    y = y_;
    z = z_;
}

void OBStar::set_velocity (double vx_, double vy_, double vz_) {
    v_x = vx_;
    v_y = vy_;
    v_z = vz_;
}

// Функция, подсчитывающия время жизни звезды на
// главной последовательности. По статье Hurley, 2000

double OBStar::get_time_on_MS () {
    return t_ms*1e6;
}

double OBStar::get_R (double T) {
    return R;
}

double OBStar::get_Z () {
    return Z;
}

double OBStar::get_L () {
    return L;
}

double OBStar::get_t_bgb() {
    return (t_bgb-t_ms)*1e6;
}

double OBStar::get_t_inf() {
    return t_inf*1e6;
}

double OBStar::get_t_HeI() {
    return (t_HeI - t_ms)*1e6;
}

double OBStar::get_M_c_SN() {
    return M_c_SN;
}

double OBStar::get_t_He() {
    return t_He*1e6;
}

double OBStar::get_c_HeI() {
    return M_c_HeI;
}

double OBStar::get_c_HeF() {
    return M_c_HeF;
}

double OBStar::get_c_BGB() {
    return M_c_BGB;
}

double OBStar::get_c_BAGB() {
    return M_c_BAGB;
}

double OBStar::get_c_DU() {
    return M_c_DU;
}
