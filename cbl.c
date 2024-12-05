/*
 * Copyright (C) 2007 Vladimir Sergeev
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  Vladimir Sergeev (vladimir.segeev [at] gmail.com) */


#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_ieee_utils.h>
#include "gsl_sum_borel.h"


/* n > k > 0 */
static int
c_n_k (int n, int k)
{
    int i, r = 1;

    for (i = 1; i <= k; i++) {
        r *= (n - k + i);
        r /= i;
    }

    return r;
}


gsl_sum_borel_leroy_workspace *
gsl_sum_borel_leroy_alloc(const size_t n)
{
    gsl_sum_borel_leroy_workspace *w;

    if (n == 0)
        GSL_ERROR_VAL ("workspace length n must be positive integer",
            GSL_EDOM, 0);

    w = (gsl_sum_borel_leroy_workspace *)
          malloc(sizeof(gsl_sum_borel_leroy_workspace));

    if (w == 0)
        GSL_ERROR_VAL("failed to allocate space for workspace struct",
            GSL_ENOMEM, 0);

    w->bseries = (double *) malloc(n * sizeof(double));

    if (w->bseries == 0) {
        free (w);

        GSL_ERROR_VAL("failed to allocate space for B(orel)L(eroy)  series",
            GSL_ENOMEM, 0);
    }

    w->useries = (double *) malloc(n * sizeof(double));

    if (w->useries == 0) {
        free(w->bseries);
        free(w);

        GSL_ERROR_VAL("failed to allocate space for C(onfrom)B(orel)L(eroy) series",
            GSL_ENOMEM, 0);
    }

    w->capacity = n;
    w->used_size = 0;
    w->a = w->b = w->Beta = w->e = 0.0;

    return w;
}

void
gsl_sum_borel_leroy_free(gsl_sum_borel_leroy_workspace *w)
{
    free(w->bseries);
    free(w);
}

/* Complex Borel-Leroy summation method of asymptotic series 's'  where
   a is 1/R (R is radius of convergence), b and beta - Borel-Leroy bare
   parameters, e is expansion parameter of 's' series.   */
int
gsl_sum_complex_borel_leroy(const double *series,
                            const size_t series_size,
                            double e, double b, double Beta, double a,
                            double epsabs, double epsrel,
                            gsl_sum_borel_leroy_workspace *w, double *result)
{
    if (series_size > w->capacity)
        GSL_ERROR ("series size exceeds available workspace", GSL_EINVAL) ;

    w->a = a;
    w->b = b;
    w->Beta = Beta;
    w->e = e;

    if (series_size == 0)
    {
        *result = 0.0;
        w->used_size = 0;
        return GSL_SUCCESS;
    }
    else if (series_size == 1)
    {
        *result = series[0];
        w->used_size = 1;
        return GSL_SUCCESS;
    }
    else
    {
        int n, m;
        double r;

        double gamma_arg = b + 1.0;
        double r0 = (4.0 / a);

        /* 0-order term */
        w->useries[0] = w->bseries[0] = series[0] / gsl_sf_gamma(gamma_arg);
        w->used_size = 1;

        /* n-order term of Borel-Leroy series, 1 < n < N */
        for (n = 1; n < series_size; n++)
        {
            gamma_arg += Beta;
            w->bseries[n] = series[n] / gsl_sf_gamma(gamma_arg);
            w->used_size++;

            /* n-order of conform Borel-Leroy series */
            r = r0;
            w->useries[n] = 0;
            for (m = 1; m <= n; m++)
            {
                w->useries[n] += w->bseries[m] * r * c_n_k(n + m - 1, n - m);
                r *= r0;
            }
        }

        /* Finally, we integrate the Borel transform:
	    \int_{0}^{\infty} e^{-t} t^{b} \sum_{k=0}^{N-1} U[k] dt
           */

        gsl_ieee_env_setup();

        double error;
        gsl_function F;
        F.function = gsl_sum_complex_borel_leroy_integrand;
        F.params = w;

        gsl_integration_workspace *iw =
		gsl_integration_workspace_alloc(1000);

        gsl_integration_qagiu(&F, 0, epsabs, epsrel, 1000, iw, result, &error);

        gsl_integration_workspace_free(iw);

    }

    return GSL_SUCCESS;
}


double
gsl_sum_complex_borel_leroy_integrand(double t, void * params)
{
    gsl_sum_borel_leroy_workspace * p =
        (gsl_sum_borel_leroy_workspace *) params;

    double a = p->a;
    double b = p->b;
    double Beta = p->Beta;
    double e = p->e;
    double * s = p->useries;
    int N = p->used_size;

    int k;

    double tmp = sqrt(1.0 + a * e * pow(t, Beta));
    double uet = (tmp - 1.0) / (tmp + 1.0);

    double o = uet;
    double r = s[0];
    for (k = 1; k < N; k++) {
        r += s[k] * o;
        o *= uet;
    }

    double f = pow(t, b) * exp(-t) * r;

    return f;
}
