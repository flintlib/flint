/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#ifndef DOUBLE_EXTRAS_H
#define DOUBLE_EXTRAS_H

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#include <float.h>
#include <gmp.h>
#include "flint.h"

#define ulong mp_limb_t

#ifdef __cplusplus
 extern "C" {
#endif

#define D_BITS 53
#define D_EPS 2.2204460492503130808e-16
#define D_INF HUGE_VAL
#define D_NAN (HUGE_VAL - HUGE_VAL)

double d_randtest(flint_rand_t state);

static __inline__ double
d_polyval(const double * poly, int len, double x)
{
    double t;
    int i;

    for (t = poly[len-1], i = len-2; i >= 0; i--)
        t = poly[i] + x * t;

    return t;
}

double d_lambertw(double x);

#ifdef __cplusplus
}
#endif

#endif
