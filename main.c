#include <stdio.h>
#include <math.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_ieee_utils.h>
#include "gsl_sum_borel.h"

#include "series.h"
#include "series.c"

#include "pade.h"
#include "pade.c"

#define EPSILON (1.0)
#define N       (1.0)
#define BETA    (1.0)
#define B0      (3.0 + 0.5 * N)
#define B       (B0 + 1.5)
#define A       (3.0/(N + 8.0))

#define ZETA_OF_3  (1.202056903)
#define H          (1.726092446232)
#define C1         (-0.493924834455)
#define C2         (-0.251106206540)
#define C3         (-0.169940489787)
#define C4         (1.806477275250)
#define A1         (-3.0/8.0)
#define A2         (-15.0/64.0)
#define A3         (-5.0/32.0)
#define A4         (45.0/32.0)

void print_series(char *prefix, double series[]);
void fill_series(double s[], double n);

void fill_series_ex(double s[], double n);

int main(int argc, char **argv)
{
    // Lev's values
    double z4[5];// = { 0.0, 0.000000000000, 0.013446156412, 0.011036395410, -0.005592766203 };
    double z3[5];// = { 0.0, 0.000000000000, 0.013446156412, 0.011036395410, 0.000000000000 };
    double z2[5];// = { 0.0, 0.000000000000, 0.013446156412, 0.000000000000, 0.000000000000 };
    double z1[5] = { 0.0, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000 };
    double z0[5] = { 0.0, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000 };

    fill_series (z4, N);
    fill_series (z3, N);
    fill_series (z2, N);
    z4[1] = z4[0] = z3[4] = z3[1] = z3[0] = z2[4] = z2[3] = z2[1] = z2[0] = 0.0;
//    fill_series (z1, N);
//    fill_series (z0, N);

    /* E-EXPANSION CM PADE */
    double EPSILON_SERIES[5];
    double CONFORM_SERIES[5];
    double PADE_4_1_A_SERIES[5];
    double PADE_3_2_A_SERIES[5];
    double PADE_2_3_A_SERIES[5];
    double PADE_3_1_SERIES[5];
    double PADE_2_2_SERIES[5];

    int i;


    gsl_sum_borel_leroy_workspace *w = gsl_sum_borel_leroy_alloc(100);

    /* Compute e-expansion of z for in all known orders L <= 4 */
    EPSILON_SERIES[0] = z4[0];
    EPSILON_SERIES[1] = z4[0] + EPSILON * z4[1];
    EPSILON_SERIES[2] = z4[0] + EPSILON * (z4[1] + EPSILON * z4[2]);
    EPSILON_SERIES[3] = z4[0] + EPSILON * (z4[1] + EPSILON * (z4[2] + EPSILON * z4[3]));
    EPSILON_SERIES[4] = z4[0] + EPSILON * (z4[1] + EPSILON * (z4[2] + EPSILON * (z4[3] + EPSILON * z4[4])));

    print_series("E-EXPANSION", EPSILON_SERIES);

    /* Conform Mapping */
    gsl_sum_complex_borel_leroy(z0, 1, EPSILON, B, BETA, A, 0.0, 1e-7, w, &CONFORM_SERIES[0]);
    gsl_sum_complex_borel_leroy(z1, 2, EPSILON, B, BETA, A, 0.0, 1e-7, w, &CONFORM_SERIES[1]);
    gsl_sum_complex_borel_leroy(z2, 3, EPSILON, B, BETA, A, 0.0, 1e-7, w, &CONFORM_SERIES[2]);
    gsl_sum_complex_borel_leroy(z3, 4, EPSILON, B, BETA, A, 0.0, 1e-7, w, &CONFORM_SERIES[3]);
    gsl_sum_complex_borel_leroy(z4, 5, EPSILON, B, BETA, A, 0.0, 1e-7, w, &CONFORM_SERIES[4]);

    print_series("CM", CONFORM_SERIES);

    for (i = 0; i < 5; ++i)
    {
        z0[i] = (i < 1) ? w->bseries[i] : 0.0;
        z1[i] = (i < 2) ? w->bseries[i] : 0.0;
        z2[i] = (i < 3) ? w->bseries[i] : 0.0;
        z3[i] = (i < 4) ? w->bseries[i] : 0.0;
        z4[i] = (i < 5) ? w->bseries[i] : 0.0;
    }

    printf("(%.12f)*(e*t)^2 + (%.12f)*(e*t)^3 + (%.12f)*(e*t)^4\n", z4[2], (z4[3]+A*z4[2]), (z4[4]+A*z4[3]));
    printf("------------------------------------------------------------------------------\n");
    printf("                       1 + (%.3f)*e*t\n", A);

    // Pade-Borel-Leroy with approx [4/1] with singularity
    gsl_sum_pade_borel_a(z2, EPSILON, B, BETA, A, 4, 1, 0.0, 1e-7, &PADE_4_1_A_SERIES[2]);
    gsl_sum_pade_borel_a(z3, EPSILON, B, BETA, A, 4, 1, 0.0, 1e-7, &PADE_4_1_A_SERIES[3]);
    gsl_sum_pade_borel_a(z4, EPSILON, B, BETA, A, 4, 1, 0.0, 1e-7, &PADE_4_1_A_SERIES[4]);

    print_series("[4,1;a]", PADE_4_1_A_SERIES);


    // Pade-Borel-Leroy with approx [3/2] with singularity
//    gsl_sum_pade_borel_a(z2, EPSILON, B, BETA, A, 3, 2, 0.0, 1e-7, &PADE_3_2_A_SERIES[2]);
//    gsl_sum_pade_borel_a(z3, EPSILON, B, BETA, A, 3, 2, 0.0, 1e-7, &PADE_3_2_A_SERIES[3]);
//    gsl_sum_pade_borel_a(z4, EPSILON, B, BETA, A, 3, 2, 0.0, 1e-7, &PADE_3_2_A_SERIES[4]);

//    print_series("[3,2;a]", PADE_3_2_A_SERIES);

    // Pade-Borel-Leroy with approx [2/3] with singularity
    gsl_sum_pade_borel_a(z2, EPSILON, B, BETA, A, 2, 3, 0.0, 1e-7, &PADE_2_3_A_SERIES[2]);
    gsl_sum_pade_borel_a(z3, EPSILON, B, BETA, A, 2, 3, 0.0, 1e-7, &PADE_2_3_A_SERIES[3]);
    gsl_sum_pade_borel_a(z4, EPSILON, B, BETA, A, 2, 3, 0.0, 1e-7, &PADE_2_3_A_SERIES[4]);

    print_series("[2,3;a]", PADE_2_3_A_SERIES);

    // Pade-Borel-Leroy with approx [2/2]
    gsl_sum_pade_borel(z2, EPSILON, B, BETA, A, 2, 2, 0.0, 1e-7, &PADE_2_2_SERIES[2]);
    gsl_sum_pade_borel(z3, EPSILON, B, BETA, A, 2, 2, 0.0, 1e-7, &PADE_2_2_SERIES[3]);
    gsl_sum_pade_borel(z4, EPSILON, B, BETA, A, 2, 2, 0.0, 1e-7, &PADE_2_2_SERIES[4]);

    print_series("[2,2]", PADE_2_2_SERIES);

    // Pade-Borel-Leroy with approx [3/1]
    gsl_sum_pade_borel(z2, EPSILON, B, BETA, A, 3, 1, 0.0, 1e-7, &PADE_3_1_SERIES[2]);
    gsl_sum_pade_borel(z3, EPSILON, B, BETA, A, 3, 1, 0.0, 1e-7, &PADE_3_1_SERIES[3]);
    gsl_sum_pade_borel(z4, EPSILON, B, BETA, A, 3, 1, 0.0, 1e-7, &PADE_3_1_SERIES[4]);

    print_series("[3,1]", PADE_2_2_SERIES);

    printf("epsilon  b        a         n\n");
    printf("%f %f %.7f %f\n", EPSILON, B, A, N);

    printf("L Eps-Exp   CM        [4/1;a]   [3/2;a]   [2/3;a]   [2/2]     [3/1]\n");
    for (i = 2; i < 5; i++)
    {

            printf("%d %.7f %.7f %.7f %.7f %.7f %.7f %.7f\n",
                i,
                EPSILON_SERIES[i],
                CONFORM_SERIES[i],
                PADE_4_1_A_SERIES[i],
                PADE_3_2_A_SERIES[i],
                PADE_2_3_A_SERIES[i],
                PADE_2_2_SERIES[i],
                PADE_3_1_SERIES[i]
               );
    }


    gsl_sum_borel_leroy_free(w);

    return 0;
}

void print_series(char * prefix, double series[])
{
    printf("%s - done.\n", prefix);

    int i;

    printf("\n");

    for (i = 0; i < 5; i++)
        printf("%s %d %.10f\n", prefix, i, series[i]);

}

void fill_series(double s[], double n)
{
    if (n == 1.0)
    {
        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 0.01344615641;
        s[3] = 0.01103639541;
        s[4] =-0.005580746254;
    }
    else if (n == 2.0)
    {
        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 0.01452184892;
        s[3] = 0.01105875302;
        s[4] =-0.005267111038;
    }
    else if (n == 3.0)
    {
        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 0.01500191004;
        s[3] = 0.01053165601;
        s[4] =-0.004977808175;
    }
    else
    {
        fill_series_ex(s, n);
    }
}

void fill_series_ex(double s[], double n)
{
    double k0 = (3.0*n + 14.0)/17.0;
    double k1 = (n + 8.0)/9.0;
    double k2 = (n*n + 6.0*n + 20.0)/27.0;
    double k3 = (n + 2.0)/3.0;
    double k4 = (5.0*n + 22.0)/27.0;

    double u1 = 2.0/(3.0*k1);
    double u2 = 34.0*k0/(81.0*k1*k1*k1);
    double u3 = 2.0*((-33.0*n*n*n + 110.0*n*n + 1760.0*n + 4544.0)/5832.0 - 4.0*ZETA_OF_3*k1*k4)/(27.0*k1*k1*k1*k1*k1);

    double r = (n + 2.0)/72.0; // k3/24

/*
    double d = (-33.0*n*n*n + 110.0*n*n + 1760.0*n + 4544.0)/5832.0 - 4.0*ZETA_OF_3*k1*k4;
    double j = 27.0 * pow (k1, 5.0);
    double jj = 1.0/j;
    double f = 2.0*k0 * jj;
*/
    // Fix: n=2 k1=1.11111
    if (n == 2.0)
        u3 = -0.2151539104;

/*
    printf("k0=%f\nk1=%f\nk2=%f\nk3=%f\nk4=%f\n", k0, k1, k2, k3, k4);
    printf("u1=%f\nu2=%f\nu3=%f\n", u1, u2, u3);
    printf("d = %f\nf = %f\nj = %f\njj = %f\n",d,f,j,jj);
*/

    s[0] = 0.0;//2.0;
    s[1] = 0.0;
    s[2] = r*u1*u1*(H - 1.0);
    s[3] = r*(2.0*u1*u2*(H - 1.0) + u1*u1*u1*(H*C1 - A1)*k1);
    s[4] = r*((2.0*u1*u3 + u2*u2)*(H - 1.0) + 3.0*u1*u1*u2*(H*C1 - A1)*k1 + u1*u1*u1*u1*((H*C2 - A2)*k2 + (H*C3 - A3)*k3 + (H*C4 - A4)*k4));
}
