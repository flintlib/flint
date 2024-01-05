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
_arb_exp_taylor_naive(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N)
{
    ulong k;
    mp_ptr s, t, u, v;
    mp_size_t nn = xn + 1;

    if (N == 0)
    {
        flint_mpn_zero(y, xn + 1);
        error[0] = 0;
        return;
    }

    s = flint_malloc(sizeof(mp_limb_t) * (nn + 1));
    t = flint_malloc(sizeof(mp_limb_t) * nn);
    v = flint_malloc(sizeof(mp_limb_t) * nn);
    u = flint_malloc(sizeof(mp_limb_t) * 2 * nn);

    /* s = 1 */
    flint_mpn_zero(s, nn);
    s[nn] = 1;

    /* t = v = x */
    flint_mpn_zero(t, nn);
    flint_mpn_copyi(t + 1, x, xn);
    flint_mpn_copyi(v, t, nn);

    for (k = 1; k < N; k++)
    {
        s[nn] += mpn_add_n(s, s, t, nn);

        /* t = t * x / (k + 1) */
        flint_mpn_mul_n(u, t, v, nn);
        flint_mpn_copyi(t, u + nn, nn);
        mpn_divrem_1(t, 0, t, nn, k + 1);
    }

    flint_mpn_copyi(y, s + 1, xn + 1);
    error[0] = 2;

    flint_free(s);
    flint_free(t);
    flint_free(u);
    flint_free(v);
}

