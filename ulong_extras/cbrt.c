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

    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#include "flint.h"
#include "ulong_extras.h"

/* function to get a good estimate of the cube root */

double
cbrt_estimate(double a)
{ 
   slong i = *((slong *) &a);
   const slong s = (1 << 10) - 1;
   i = (i - (s << 52)) * 0.333333333333333 + (s << 52); /* 0.333333333333333 = 1/3*/
   return *((double *) &i);
}

/* failsafe cube root function using pow and checking for overflow */

mp_limb_t
n_cbrt_failsafe(mp_limb_t n)
{
    mp_limb_t val, upper_limit;
    upper_limit = 1626 + 1;     /* 1626 < (2^32)^(1/3) */

#if FLINT64
    upper_limit = 2642245 + 1;  /* 2642245 < (2^64)^(1/3) */
#endif

    val = (mp_limb_t) pow((double)n, 1./3.) + 1;

    if (val < upper_limit)
        val -= (val*val*val > n);
    else
        val -= 1;

    return val;
}

/* fast cube root function. If value returned calculated after cbrt_estimate
    and one round of newton iteration causes overflow, n_cbrt_failsafe is called */

mp_limb_t
n_cbrt(mp_limb_t n)
{
    mp_limb_t ret, upper_limit;
    double val, x, dx, xsq;
    int iter;
    iter = 2;
    val = (double)n;

    upper_limit = 1626 + 1;     /* 1626 < (2^32)^(1/3) */
#if FLINT64
    upper_limit = 2642245 + 1;  /* 2642245 < (2^64)^(1/3) */
#endif

    x = cbrt_estimate((double)n);

    /* newton iteration */

    while (iter--)
    {
        xsq = x*x;
        dx = val / xsq;
        dx -= x;
        dx*= 0.333333333333333;     /* dx = dx * (1/3) */     
        x+=dx;
    }
    ret = x;
    if (ret>= upper_limit)      
    {
        ret = n_cbrt_failsafe(n);
        return ret;
    }

    while (ret*ret*ret <= n)
        (ret) += 1;
    while (ret*ret*ret > n)
        (ret) -= 1;

    return ret;
}