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

#include "ulong_extras-impl.h"

mp_limb_t
n_root(mp_limb_t n, mp_limb_t root)
{
    mp_limb_t x, currval, base, upper_limit; 
    double dx;

    if (!n || !root)
        return 0;

    if (root == 1)
        return n;

    if (root == 2)
        return n_sqrt(n);

    if (root == 3)
        return n_cbrt(n);

    if (root >= FLINT_BITS || (UWORD(1) << root) > n)
        return 1;

    upper_limit = max_base[root - 1];      /* n <= upper_limit^root */

    x = n_root_estimate((double)n, root);

    /* one round of newton iteration */

    currval = n_pow(x, root-1);
    dx = n / currval;
    dx -= x;
    dx *= inv_table[root - 1];
    dx = floor(dx);
    x += dx;

    base = x;

    if (base >= upper_limit)
        base = upper_limit - 1;

    currval = n_pow(base, root);
    if (currval == n)
        goto final;

    while (currval <= n)
    {
        (base) += 1;
        currval = n_pow(base, root);
        if (base == upper_limit)
            break;
    }

    while (currval > n)
    {
        (base) -= 1;
        currval = n_pow(base, root);
    }

    final:
    return base;
}
