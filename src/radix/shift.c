/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/* Note: both lshift and rshift currently iterate in the convenient direction
   for the algorithm. The reverse direction (used by their mpn counterparts)
   would be better in terms of allowing in-place shifts with a multi-limb
   offset. */

ulong
radix_lshift_digits(nn_ptr res, nn_srcptr a, slong n, unsigned int e, const radix_t radix)
{
    FLINT_ASSERT(e >= 1 && e < radix->exp);

    ulong cy, hi, lo;
    slong i;

    cy = 0;

    ulong be = radix->bpow[e];
    ulong be2 = radix->bpow[radix->exp - e];

    n_div_precomp_t pre;
    *pre = radix->bpow_div[radix->exp - e];

    if (pre->m == 0)
    {
        for (i = 0; i < n; i++)
        {
            hi = n_divrem_precomp_m0(&lo, a[i], be2, pre);
            res[i] = cy + lo * be;
            cy = hi;
        }
    }
    else if (pre->c == 0)
    {
        for (i = 0; i < n; i++)
        {
            hi = n_divrem_precomp_c0(&lo, a[i], be2, pre);
            res[i] = cy + lo * be;
            cy = hi;
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            hi = n_divrem_precomp_c1_bounded(&lo, a[i], be2, pre);
            res[i] = cy + lo * be;
            cy = hi;
        }
    }

    return cy;
}

ulong
radix_rshift_digits(nn_ptr res, nn_srcptr a, slong n, unsigned int e, const radix_t radix)
{
    FLINT_ASSERT(e >= 1 && e < radix->exp);

    ulong cy, hi, lo;
    slong i;

    cy = 0;

    ulong be = radix->bpow[e];
    ulong be2 = radix->bpow[radix->exp - e];

    n_div_precomp_t pre;
    *pre = radix->bpow_div[e];

    if (pre->m == 0)
    {
        for (i = n - 1; i >= 0; i--)
        {
            hi = n_divrem_precomp_m0(&lo, a[i], be, pre);
            res[i] = hi + cy * be2;
            cy = lo;
        }
    }
    else if (pre->c == 0)
    {
        for (i = n - 1; i >= 0; i--)
        {
            hi = n_divrem_precomp_c0(&lo, a[i], be, pre);
            res[i] = hi + cy * be2;
            cy = lo;
        }
    }
    else
    {
        for (i = n - 1; i >= 0; i--)
        {
            hi = n_divrem_precomp_c1_bounded(&lo, a[i], be, pre);
            res[i] = hi + cy * be2;
            cy = lo;
        }
    }

    return cy * be2;
}

