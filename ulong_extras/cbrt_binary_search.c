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
n_cbrt_binary_search(mp_limb_t x)
{
    mp_limb_t low, high, mid, p, upper_limit;

    /* upper_limit is the max cube root possible for one word */

#ifdef FLINT64
    upper_limit = 2642245;  /* 2642245 < (2^64)^(1/3) */
#else
    upper_limit = 1625;     /* 1625 < (2^32)^(1/3) */
#endif

    low = 0;
    high = UWORD(1) << ((FLINT_BIT_COUNT(x) + 2) / 3);

    if (high > upper_limit) /* cube cannot be greater than upper_limit */
        high = upper_limit;

    /* binary search for cube root */
        
    while (low < high)
    {
        mid = (high + low) / 2;
        p = mid + 1;
        p = p * p * p;
        if (p == x)
            return mid + 1;
        else if (p > x)
            high = mid;
        else
            low = mid + 1;
    }
    return low;
}
