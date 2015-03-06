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
    Copyright (C) 2015 Fredrik Johansson
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
    and two iterations causes overflow, n_cbrt_failsafe is called */

mp_limb_t
n_cbrt(mp_limb_t n)
{
    int iter, bits;
    double val, x, xcub, num, den;
    mp_limb_t ret, upper_limit, low, high, mid, p;

    val = (double)n;
    bits = FLINT_BIT_COUNT(n); 
    upper_limit = 1626;     /* 1626 < (2^32)^(1/3) */
#if FLINT64
    upper_limit = 2642245;  /* 2642245 < (2^64)^(1/3) */
#endif

    if (bits < 8)           /* binary search in case of small n */
    {
        low = 0;
        high = UWORD(1) << ((FLINT_BIT_COUNT(n) + 2) / 3);
        while (low < high)
        {
            mid = low + (high - low) / 2;
            p = mid + 1;
            p = p*p*p;

            if (p == n)
                return mid + 1;
            else if (p > n)
                high = mid;
            else
                low = mid + 1;
        }
        return low;
    }
    if (bits < 46)      /* one iteration seems to be sufficient for n < 2^46 */
        iter = 1;
    else
        iter = 2;       /* 2 gives us a precise enough answer for any mp_limb_t */
    
    x = cbrt_estimate((double)n);   /* initial estimate */
   
    /* Kahan's iterations to get cube root */

    val = (double)n;
    while(iter--)
    {
        xcub = x*x*x;
        num = (xcub - val)*x;
        den = (xcub + xcub + val);
        x -= (num/den);
    }

    ret = x;
    if (ret>= upper_limit)      
    {
        if (ret == upper_limit)
            return ret;
        return n_cbrt_failsafe(n);
    }

    while (ret*ret*ret <= n)
        (ret) += 1;
    while (ret*ret*ret > n)
        (ret) -= 1;

    return ret;
}
