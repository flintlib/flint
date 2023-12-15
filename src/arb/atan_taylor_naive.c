/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mpn_extras.h"

void
_arb_atan_taylor_naive(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N, int alternating)
{
    ulong k;
    mp_ptr s, t, x1, x2, u;
    mp_size_t nn = xn + 1;

    if (N == 0)
    {
        flint_mpn_zero(y, xn);
        error[0] = 0;
        return;
    }

    if (N == 1)
    {
        flint_mpn_copyi(y, x, xn);
        error[0] = 0;
    }

    s = flint_malloc(sizeof(mp_limb_t) * nn);
    t = flint_malloc(sizeof(mp_limb_t) * nn);
    u = flint_malloc(sizeof(mp_limb_t) * 2 * nn);
    x1 = flint_malloc(sizeof(mp_limb_t) * nn);
    x2 = flint_malloc(sizeof(mp_limb_t) * nn);

    flint_mpn_zero(s, nn);
    flint_mpn_zero(t, nn);
    flint_mpn_zero(u, 2 * nn);
    flint_mpn_zero(x1, nn);
    flint_mpn_zero(x2, nn);

    /* x1 = x */
    flint_mpn_copyi(x1 + 1, x, xn);

    /* x2 = x * x */
    flint_mpn_mul_n(u, x1, x1, nn);
    flint_mpn_copyi(x2, u + nn, nn);

    /* s = t = x */
    flint_mpn_copyi(s, x1, nn);
    flint_mpn_copyi(t, x1, nn);

    for (k = 1; k < N; k++)
    {
        /* t = t * x2 */
        flint_mpn_mul_n(u, t, x2, nn);
        flint_mpn_copyi(t, u + nn, nn);

        /* u = t / (2k+1) */
        mpn_divrem_1(u, 0, t, nn, 2 * k + 1);

        if (alternating & k)
            mpn_sub_n(s, s, u, nn);
        else
            mpn_add_n(s, s, u, nn);
    }

    flint_mpn_copyi(y, s + 1, xn);
    error[0] = 2;

    flint_free(s);
    flint_free(t);
    flint_free(u);
    flint_free(x1);
    flint_free(x2);
}

