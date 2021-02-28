/*
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

/* return (x|y)*(-1)^tstbit(r,1) */
int _n_jacobi_unsigned(mp_limb_t x, mp_limb_t y, unsigned int r)
{
    mp_limb_t t, st;
    int e;

    FLINT_ASSERT(y & 1);

    r ^= 2;

    while (y > 1)
    {
        if (x == 0)
            return 0;

        /* x = odd part of x */
        count_trailing_zeros(e, x);
        x >>= e;
        r ^= ((y ^ (y>>1)) & (2*e)); /* (2|y) = (-1)^((y^2-1)/8) */

        /* (x, y) = (|x - y|, min(x, y)) */
        sub_ddmmss(st, t, UWORD(0), x, UWORD(0), y);
        r ^= (x & y & st);  /* if y > x, (x|y) = (y|x)*(-1)^((x-1)(y-1)/4) */
        y += (st & t);
        x = (t ^ st) - st;
    }

    return (int)(r & 2) - 1;
}

int n_jacobi_unsigned(mp_limb_t x, mp_limb_t y)
{
    return _n_jacobi_unsigned(x, y, 0);
}

int n_jacobi(mp_limb_signed_t x, mp_limb_t y)
{
    return _n_jacobi_unsigned(FLINT_ABS(x), y, FLINT_SIGN_EXT(x) & y);
}

