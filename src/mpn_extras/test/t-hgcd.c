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

/* uses _flint_mpn_gcd_test_pair and _flint_mpn_gcd_test_len (t-gcd.c) */

static void
_mpz_set_mpn(mpz_t x, mp_srcptr d, mp_size_t n)
{
    MPN_NORM(d, n);
    mpz_import(x, n, -1, sizeof(mp_limb_t), 0, 0, d);
}

/* Checks the defining properties of the half-gcd: with M the returned
   matrix and (a', b') the reduced operands, (a; b) = M (a'; b'),
   det M = 1, and a' and b' have the returned size; the reduction reaches
   (almost) half the size, and the entries of M fit in about half. */
TEST_FUNCTION_START(flint_mpn_hgcd, state)
{
    slong iter;
    mpz_t a, b, a2, b2, m[2][2], t, u;
    int i, j;

    mpz_init(a);
    mpz_init(b);
    mpz_init(a2);
    mpz_init(b2);
    mpz_init(t);
    mpz_init(u);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            mpz_init(m[i][j]);

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong n, nn;
        mp_ptr ap, bp, mp, tp;
        flint_mpn_hgcd_matrix_struct M;

        n = _flint_mpn_gcd_test_len(state) + 2;
        _flint_mpn_gcd_test_pair(a, b, state, n, n_randint(state, 4) ? n : 1 + n_randint(state, n));

        n = mpz_size(a);
        if (n < 3)
            continue;

        ap = flint_calloc(n + 1, sizeof(mp_limb_t));
        bp = flint_calloc(n + 1, sizeof(mp_limb_t));
        mp = flint_malloc(FLINT_MPN_HGCD_MATRIX_INIT_ITCH(n) * sizeof(mp_limb_t));
        tp = flint_malloc(_flint_mpn_hgcd_itch(n) * sizeof(mp_limb_t));

        flint_mpn_copyi(ap, a->_mp_d, n);
        flint_mpn_copyi(bp, b->_mp_d, mpz_size(b));
        if (n_randint(state, 2))
        {
            FLINT_SWAP(mp_ptr, ap, bp);
            mpz_swap(a, b);
        }

        _flint_mpn_hgcd_matrix_init(&M, n, mp);
        nn = _flint_mpn_hgcd(ap, bp, n, &M, tp);

        if (nn != 0)
        {
            _mpz_set_mpn(a2, ap, nn);
            _mpz_set_mpn(b2, bp, nn);

            if (nn > n || (ap[nn - 1] | bp[nn - 1]) == 0)
                TEST_FUNCTION_FAIL("bad size: n = %wd, nn = %wd\n", n, nn);

            /* reduced to |a' - b'| < B^(n/2 + 1) */
            mpz_sub(t, a2, b2);
            if (mpz_sizeinbase(t, 2) > (ulong) FLINT_BITS * (n / 2 + 1))
                TEST_FUNCTION_FAIL("not reduced: n = %wd, nn = %wd\na = %{mpz}\nb = %{mpz}\n",
                    n, nn, a, b);

            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                {
                    if (M.n > (n + 1) / 2)
                        TEST_FUNCTION_FAIL("matrix too large: n = %wd, M.n = %wd\n", n, M.n);
                    _mpz_set_mpn(m[i][j], M.p[i][j], M.n);
                }

            /* (a; b) = M (a'; b') */
            mpz_mul(t, m[0][0], a2);
            mpz_addmul(t, m[0][1], b2);
            mpz_mul(u, m[1][0], a2);
            mpz_addmul(u, m[1][1], b2);
            if (mpz_cmp(t, a) != 0 || mpz_cmp(u, b) != 0)
                TEST_FUNCTION_FAIL("M (a'; b') != (a; b): n = %wd, nn = %wd\na = %{mpz}\nb = %{mpz}\n",
                    n, nn, a, b);

            /* det M = 1 */
            mpz_mul(t, m[0][0], m[1][1]);
            mpz_submul(t, m[0][1], m[1][0]);
            if (mpz_cmp_ui(t, 1) != 0)
                TEST_FUNCTION_FAIL("det M != 1: n = %wd\n", n);
        }

        flint_free(ap);
        flint_free(bp);
        flint_free(mp);
        flint_free(tp);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(a2);
    mpz_clear(b2);
    mpz_clear(t);
    mpz_clear(u);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            mpz_clear(m[i][j]);

    TEST_FUNCTION_END(state);
}
