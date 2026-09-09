/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_bdiv_qr, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, b, q, r, c, ac;
        mp_size_t an, bn, n, cn, i;
        int alg;

        an = 1 + n_randint(state, 40);
        bn = 1 + n_randint(state, 40);
        n = 1 + n_randint(state, 60);
        if (n_randint(state, 50) == 0)
        {
            an = 1 + n_randint(state, 1200);
            bn = 1 + n_randint(state, 1200);
            n = 1 + n_randint(state, 1200);
        }

        a = flint_malloc(an * sizeof(mp_limb_t));
        b = flint_malloc(bn * sizeof(mp_limb_t));
        q = flint_malloc(n * sizeof(mp_limb_t));
        r = flint_malloc(bn * sizeof(mp_limb_t));
        cn = n + bn;
        c = flint_malloc(cn * sizeof(mp_limb_t));
        ac = flint_malloc(cn * sizeof(mp_limb_t));

        flint_mpn_rrandom(a, state, an);
        flint_mpn_rrandom(b, state, bn);
        b[0] |= 1;

        alg = n_randint(state, 4);
        if (alg == 0)
            flint_mpn_bdiv_qr(q, r, a, an, b, bn, n);
        else if (alg == 1)
            flint_mpn_bdiv_qr_classical(q, r, a, an, b, bn, n);
        else if (alg == 2)
            flint_mpn_bdiv_qr_karp_markstein(q, r, a, an, b, bn, n);
        else
        {
            flint_mpn_bdiv_qr(q, NULL, a, an, b, bn, n);
            /* the q-only version must agree with the full one */
            {
                mp_ptr q2 = flint_malloc(n * sizeof(mp_limb_t));
                flint_mpn_bdiv_qr_karp_markstein(q2, r, a, an, b, bn, n);
                if (mpn_cmp(q, q2, n) != 0)
                    TEST_FUNCTION_FAIL("q-only mismatch: an = %wd, bn = %wd, n = %wd\n", an, bn, n);
                flint_free(q2);
            }
        }

        /* check a == q b + B^n r  (mod B^(n + bn)) */
        if (n >= bn)
            flint_mpn_mul(c, q, n, b, bn);
        else
            flint_mpn_mul(c, b, bn, q, n);
        mpn_add_n(c + n, c + n, r, bn);

        for (i = 0; i < cn; i++)
            ac[i] = (i < an) ? a[i] : 0;

        if (mpn_cmp(c, ac, cn) != 0)
            TEST_FUNCTION_FAIL("alg = %d, an = %wd, bn = %wd, n = %wd\na = %{ulong*}\nb = %{ulong*}\nq = %{ulong*}\nr = %{ulong*}\n",
                alg, an, bn, n, a, an, b, bn, q, n, r, bn);

        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(r);
        flint_free(c);
        flint_free(ac);
    }

    TEST_FUNCTION_END(state);
}
