typedef struct
{
    double z0, z1, z2, z3, z4;
    double b, Beta, a, e;
} pade_borel_f_params;


int
gsl_sum_pade_borel(double *series, double e, double b, double Beta, double a,
                   int nom_degree, int denom_degree,
                   double epsabs, double epsrel, double *result);

int
gsl_sum_pade_borel_a(double *series, double e, double b, double Beta, double a,
                     int nom_degree, int denom_degree,
                     double epsabs, double epsrel, double *result);

//Pade with known singularity type
double pade_borel_4_1_a(double t, void *params);

double pade_borel_3_2_a(double t, void *params);

double pade_borel_2_3_a(double t, void *params);

/*
double pade_borel_3_1_a(double t, void *params);

double pade_borel_2_1_a(double t, void *params);

double pade_borel_1_1_a(double t, void *params);

double pade_borel_0_1_a(double t, void *params);

double pade_borel_0_4(double t, void *params);

double pade_borel_1_3(double t, void *params);
*/

double pade_borel_2_2(double t, void *params);

double pade_borel_3_1(double t, void *params);

/*
double pade_borel_0_3(double t, void *params);

double pade_borel_1_2(double t, void *params);

double pade_borel_2_1(double t, void *params);

double pade_borel_0_2(double t, void *params);

double pade_borel_1_1(double t, void *params);

double pade_borel_0_1(double t, void *params);
*/
