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
n_cbrt(mp_limb_t n)
{
    int bits;
    double val, x, xcub, num, den;
    mp_limb_t ret, upper_limit;

    /* Taking care of smaller roots */

    if (n < 125)
        return (n >= 1) + (n >= 8) + (n >= 27) + (n >= 64);
    if (n < 1331)
        return 5 + (n >= 216) + (n >= 343) + (n >= 512) + (n >= 729) + (n >= 1000);
    if (n < 4913)
        return 11 + (n >= 1728) + (n >= 2197) + (n >= 2744) + (n >= 3375) + (n >= 4096);

    val = (double)n;
    bits = FLINT_BIT_COUNT(n); 

    if (bits > 46)    /* for larger numbers chebyshev approximation method is faster */
        return n_cbrt_chebyshev_approx(n);

    /* upper_limit is the max cube root possible for one word */

#ifdef FLINT64
    upper_limit = 2642245;  /* 2642245 < (2^64)^(1/3) */
#else
    upper_limit = 1625;     /* 1625 < (2^32)^(1/3) */
#endif

    x = n_cbrt_estimate((double)n);   /* initial estimate */

    /* Kahan's iterations to get cube root */

    xcub = x * x * x;
    num = (xcub - val) * x;
    den = (xcub + xcub + val);
    x -= (num / den);
    ret = x;

    /* In case ret^3 or (ret+1)^3 will cause overflow */

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

