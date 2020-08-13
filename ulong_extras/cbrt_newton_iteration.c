/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t
n_cbrt_newton_iteration(mp_limb_t n)
{
    int iter, bits;
    mp_limb_t ret;
    double val, x, xsq, dx;

    /* upper_limit is the max cube root possible for one word */

#ifdef FLINT64
    const mp_limb_t upper_limit = 2642245;  /* 2642245 < (2^64)^(1/3) */
#else
    const mp_limb_t  upper_limit = 1625;    /* 1625 < (2^32)^(1/3) */
#endif
    
    val = (double)n;
    bits = FLINT_BIT_COUNT(n); 
    if (bits < 46)      /* one iteration seems to be sufficient for n < 2^46 */
        iter = 1;
    else
        iter = 2;       /* 2 gives us a precise enough answer for any mp_limb_t */
    
    x = n_cbrt_estimate((double)n);         /* initial estimate */

    /* Newton's iterations to get cube root */
    val = (double)n;
    while(iter--)
    {
        xsq = x * x;
        dx = val / xsq;
        dx -= x;
        dx *= 0.333333333333333;     /* dx = dx * (1/3) */     
        x += dx;
    }
    /* In case ret^3 or (ret+1)^3 will cause overflow */

    ret = x;
    if (ret >= upper_limit)      
    {
        if (n >= upper_limit * upper_limit * upper_limit)
            return upper_limit;
        ret = upper_limit - 1;
    }
    while (ret * ret * ret <= n)
    {
        (ret) += 1;
        if (ret == upper_limit)
            break;
    }
    while (ret * ret * ret > n)
        (ret) -= 1;
    return ret;
}
