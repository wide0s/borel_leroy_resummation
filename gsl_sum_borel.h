/* sum/gsl_sum_borel.h
 *
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


#ifndef __GSL_SUM_BOREL_H__
#define __GSL_SUM_BOREL_H__

#include <stdlib.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

__BEGIN_DECLS

/*  Workspace for conformal Borel-Leroy transform
 *
 *   capacity    = number of terms the workspace can handle
 *   used_size   = actual number of terms in workspace
 *   a           = 1/a is radius of convergance of base series
 *   b           = bare parameter
 *   Beta        = bare parameter
 *   e           = expansion parameter of base series
 *   bseries     = Borel-Leroy series
 *   useries     = Conformal Borel-Leroy series
 */

typedef struct
{
    size_t capacity;
    size_t used_size;
    double a;
    double b;
    double Beta;
    double e;
    double * bseries;
    double * useries;
} gsl_sum_borel_leroy_workspace;

gsl_sum_borel_leroy_workspace *gsl_sum_borel_leroy_alloc(const size_t n);

void
gsl_sum_borel_leroy_free(gsl_sum_borel_leroy_workspace *w);

/* Conformal Borel-Leroy summation method
 *
 *   series      = base asymptotic series
 *   series_size = size of series
 *   e           = base series expansion parameter
 *   b, Beta     = bare parameters
 *   a           = 1/a is radius of convergance of base series
 *   epsabs      = abs error of qagiu integration
 *   epsrel      = rel error of qagiu integration
 *   w           = workspace
 *   result      = summation result
 */

int
gsl_sum_complex_borel_leroy(const double *series,
                            const size_t series_size,
                            double e, double b, double Beta, double a,
                            double epsabs, double epsrel,
                            gsl_sum_borel_leroy_workspace *w, double *result);

double
gsl_sum_complex_borel_leroy_integrand(double t, void *params);

__END_DECLS

#endif /* __GSL_SUM_BOREL_H__ */
