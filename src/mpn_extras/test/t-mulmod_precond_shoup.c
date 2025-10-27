/*
    Copyright (C) 2013 William Hart
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_mulmod_precond_shoup, state)
{
    int i;
    slong n;
    ulong norm;

    mp_ptr a, apre, b, d, dnormed, dinv, r1, r2, t, u;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        n = 1 + n_randint(state, 10);
        norm = n_randint(state, 2);

        a = flint_malloc(sizeof(mp_limb_t) * n);
        apre = flint_malloc(sizeof(mp_limb_t) * n);
        b = flint_malloc(sizeof(mp_limb_t) * n);
        d = flint_malloc(sizeof(mp_limb_t) * n);
        dnormed = flint_malloc(sizeof(mp_limb_t) * n);
        dinv = flint_malloc(sizeof(mp_limb_t) * n);
        t = flint_malloc(sizeof(mp_limb_t) * 2 * n);
        u = flint_malloc(sizeof(mp_limb_t) * (n + 1));
        r1 = flint_malloc(sizeof(mp_limb_t) * n);
        r2 = flint_malloc(sizeof(mp_limb_t) * n);

        flint_mpn_rrandom(d, state, n);
        do
        {
            flint_mpn_rrandom(d + n - 1, state, 1);
            if (norm == 0)
                d[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));
            else
                d[n - 1] >>= 1;
        }
        while (d[n - 1] == 0);

        norm = flint_clz(d[n - 1]);
        if (norm == 0)
            mpn_copyi(dnormed, d, n);
        else
            mpn_lshift(dnormed, d, n, norm);
        flint_mpn_preinvn(dinv, dnormed, n);

        /* reduce a, b mod d */
        flint_mpn_rrandom(a, state, n);
        mpn_tdiv_qr(t, a, 0, a, n, d, n);
        flint_mpn_rrandom(b, state, n);
        mpn_tdiv_qr(t, b, 0, b, n, d, n);

        mpn_mul_n(t, a, b, n);
        mpn_tdiv_qr(u, r1, 0, t, 2 * n, d, n);

        flint_mpn_mulmod_precond_shoup_precompute(apre, a, n, dnormed, dinv, norm);
        flint_mpn_mulmod_precond_shoup(r2, a, apre, b, n, d, norm);

        if (mpn_cmp(r1, r2, n))
        {
            flint_printf("FAIL\n");
            flint_printf("n = %wd, norm = %wu\n", n, norm);
            flint_printf("d = "); flint_mpn_debug(d, n);
            flint_printf("a = "); flint_mpn_debug(a, n);
            flint_printf("apre = "); flint_mpn_debug(apre, n + (norm == 0));
            flint_printf("b = "); flint_mpn_debug(b, n);
            flint_printf("r1 = "); flint_mpn_debug(r1, n);
            flint_printf("r2 = "); flint_mpn_debug(r2, n);
            flint_abort();
        }

        flint_free(a);
        flint_free(apre);
        flint_free(b);
        flint_free(d);
        flint_free(dnormed);
        flint_free(dinv);
        flint_free(t);
        flint_free(u);
        flint_free(r1);
        flint_free(r2);
    }

    TEST_FUNCTION_END(state);
}
