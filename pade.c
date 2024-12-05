int
gsl_sum_pade_borel_a(double * series, double e, double b, double Beta, double a,
                     int nom_degree, int denom_degree,
                     double epsabs, double epsrel, double *result)
{
    double error;
    gsl_function F;
    pade_borel_f_params p;


    if (((nom_degree + denom_degree) > 5) || (nom_degree < 0) || (denom_degree < 0))
    {
        GSL_ERROR("inavlid pade selected", GSL_EINVAL) ;
    }

    p.z0 = series[0];
    p.z1 = series[1];
    p.z2 = series[2];
    p.z3 = series[3];
    p.z4 = series[4],
    p.b = b, p.Beta = Beta, p.a = a, p.e = e;

    F.params = &p;

    // Pade[4,1]
    if ((nom_degree == 4) && (denom_degree == 1))
    {
        F.function = pade_borel_4_1_a;
    }
    else if ((nom_degree == 3) && (denom_degree == 2))
    {
       F.function = pade_borel_3_2_a;
    }
    else if ((nom_degree == 2) && (denom_degree == 3))
    {
       F.function = pade_borel_2_3_a;
    }
    else
    {
        GSL_ERROR("inavlid pade selected", GSL_EINVAL);
    }

    gsl_ieee_env_setup();

    gsl_integration_workspace *w =
         gsl_integration_workspace_alloc(1000);

    gsl_integration_qagiu(&F, 0, epsabs, epsrel, 1000, w, result, &error);

    gsl_integration_workspace_free(w);

    return GSL_SUCCESS;
}

int
gsl_sum_pade_borel(double *series, double e, double b, double Beta, double a,
                   int nom_degree, int denom_degree,
                   double epsabs, double epsrel, double *result)
{
    double error;
    gsl_function F;
    pade_borel_f_params p;


    if (((nom_degree + denom_degree) > 4) || (nom_degree < 0) || (denom_degree < 0))
    {
        GSL_ERROR("inavlid pade selected", GSL_EINVAL) ;
    }

    p.z0 = series[0];
    p.z1 = series[1];
    p.z2 = series[2];
    p.z3 = series[3];
    p.z4 = series[4],
    p.b = b, p.Beta = Beta, p.a = a, p.e = e;

    F.params = &p;

    if ((nom_degree == 2) && (denom_degree == 2))
    {
        F.function = pade_borel_2_2;
    }
    else if ((nom_degree == 3) && (denom_degree == 1))
    {
        F.function = pade_borel_3_1;
    }
    else
    {
        GSL_ERROR("inavlid pade selected", GSL_EINVAL);
    }

    gsl_ieee_env_setup();

    gsl_integration_workspace *w =
         gsl_integration_workspace_alloc(1000);

    gsl_integration_qags(&F, 0, 100, epsabs, epsrel, 1000, w, result, &error);

    gsl_integration_workspace_free(w);

    return GSL_SUCCESS;
}

double pade_borel_4_1_a(double t, void *params)
{
    pade_borel_f_params * p =
        (pade_borel_f_params *) params;

    double z2 = p->z2, z3 = p->z3, z4 = p->z4;
    double e = p->e, b = p->b, Beta = p->Beta, a = p->a;

    double et = e * pow(t, Beta); // e*t^beta

    double nominator, denominator, F;

    nominator = (((z4 + a*z3)*et + z3 + a*z2)*et + z2)*et*et;

    denominator = 1.0 + a*et;

    F = pow(t, b) * exp(-t) * (nominator / denominator);

    return F;
}

double pade_borel_3_2_a(double t, void * params)
{
    pade_borel_f_params * p =
        (pade_borel_f_params *) params;

    double z2 = p->z2, z3 = p->z3, z4 = p->z4;
    double e = p->e, b = p->b, Beta = p->Beta, a = p->a;
    double s;
    double et = e * pow(t, Beta); // e*t^beta

    double nominator, denominator, F;

    if ((z4 == 0.0) && (z3 == 0.0))
    {
        nominator = z2*et*et;
        denominator = 1.0;
    }
    else
    {
       s = a*z2+z3;

       nominator = pow(et, 2) * (z2 + (a*a*z2*z3-z2*z4+z2*z2)*et/s);

       //nominator = (z2*(a*z2 + z3) + (a*z2*(a*z2 + z3) + z3*z3 - z2*z4)*et)*et*et;
       //nominator = (a*z2*z2+z2*z3)*et*et + (a*z2*z3-z2*z4+z2*z2)*et*et*et;
       denominator = (1.0 + a*et) * (1.0 - (a*z3 + z4)*et/s);
    }

    F = pow(t, b) * exp(-t) * (nominator / denominator);

    return F;
}


double
pade_borel_2_3_a (double t, void * params)
{
    pade_borel_f_params * p =
        (pade_borel_f_params *) params;

    double z2 = p->z2, z3 = p->z3, z4 = p->z4;
    double e = p->e, b = p->b, Beta = p->Beta, a = p->a;

    double et = e * pow(t, Beta); // e*t^beta

    double nominator, denominator, F;

    nominator = z2*z2*z2*et*et;

    denominator = (1.0 + a*et) * ((((a*z2*(a*z2 + z3) + z3*z3 - z2*z4)*et - (a*z2*z2 + z3*z2))*et) + z2*z2);

//    denominator = (1.0 + a*et)*(z2*z2 - (a*z2+z3)*et + (a*a*z2*z2 + a*z2*z3 + z3*z3 - z4*z2)*et*et);

    F = pow(t, b) * exp(-t) * (nominator / denominator);

    return F;
}



double pade_borel_2_2(double t, void * params)
{
   pade_borel_f_params * p =
        (pade_borel_f_params *) params;

    double z0 = p->z0, z1 = p->z1, z2 = p->z2, z3 = p->z3, z4 = p->z4;
    double e = p->e, b = p->b, Beta = p->Beta;

    double et = e*pow(t, Beta); // e*t^beta

    double z1p2 = pow(z1, 2); // z1^2
    double z2p2 = pow(z2, 2); // z2^2
    double z3p2 = pow(z3, 2); // z3^2

    double z2p3 = pow(z2, 3); // z2^3

    double nominator, denominator, F;

    if ((z4 == 0.0) && (z3 == 0.0) && (z2 == 0.0) && (z1 == 0.0))
    {
        nominator = z0;
        denominator = 1.0;
    }
    else
    {
        nominator = (2.0*z2*z1*z3 - z2p3 - z1p2*z4 - z0*z3p2 + z0*z2*z4)*pow (et, 2) +
                    (z1p2*z3 - z1*z2p2 - z0*z1*z4 + z0*z2*z3)*et +
                    z0*z1*z3 - z0*z2p2;

        denominator = (-z3p2 + z2*z4)*pow (et, 2) +
                      (-z1*z4 + z2*z3)*et +
                      z1*z3 - z2p2;
    }

    F = pow(t, b) * exp(-t) * (nominator / denominator);

    return F;
}

double pade_borel_3_1(double t, void * params)
{
   pade_borel_f_params * p =
        (pade_borel_f_params *) params;

    double z0 = p->z0, z1 = p->z1, z2 = p->z2, z3 = p->z3, z4 = p->z4;
    double e = p->e, b = p->b, Beta = p->Beta;

    double et = e * pow(t, Beta); // e*t^beta

    double z3p2 = pow(z3, 2); // z3^2

    double nominator, denominator, F;

    if (z4 == 0.0)
    {
        nominator = z3 * pow(et, 3) + z2 * pow(et, 2) + z1 * et + z0;

        denominator = 1.0;
    }
    else
    {
        nominator = (z3p2 - z2*z4)*pow (et, 3) + (-z1*z4 + z2*z3)*pow (et, 2) +
                    (z1*z3 - z4*z0)*et + z3*z0;

        denominator = z3 - z4*et;
    }

    F = pow(t, b) * exp(-t) * (nominator / denominator);

    return F;
}
