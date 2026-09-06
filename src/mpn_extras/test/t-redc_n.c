/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "test_helpers.h"
#include "mpn_extras.h"

/* checks t == X * B^(-n) mod m and t < m, via t * B^n ≡ X (mod m) */
static int
redc_ok(nn_srcptr t, nn_srcptr X, nn_srcptr m, mp_size_t n)
{
    mpz_t tz, Xz, mz;
    int ok;

    if (mpn_cmp(t, m, n) >= 0)
        return 0;

    mpz_init(tz);
    mpz_init(Xz);
    mpz_init(mz);
    mpz_import(tz, n, -1, sizeof(ulong), 0, 0, t);
    mpz_import(Xz, 2 * n, -1, sizeof(ulong), 0, 0, X);
    mpz_import(mz, n, -1, sizeof(ulong), 0, 0, m);

    mpz_mul_2exp(tz, tz, n * FLINT_BITS);
    mpz_sub(tz, tz, Xz);
    ok = mpz_divisible_p(tz, mz);

    mpz_clear(tz);
    mpz_clear(Xz);
    mpz_clear(mz);
    return ok;
}

TEST_FUNCTION_START(flint_mpn_redc_n, state)
{
    for (slong iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        mp_size_t n = 1 + n_randint(state, n_randint(state, 10) == 0 ?
                                           400 : 60);
        nn_ptr m, minv, v, X, t, scratch;
        slong j;

        m = FLINT_ARRAY_ALLOC(n, ulong);
        minv = FLINT_ARRAY_ALLOC(n, ulong);
        v = FLINT_ARRAY_ALLOC(n, ulong);
        X = FLINT_ARRAY_ALLOC(2 * n, ulong);
        t = FLINT_ARRAY_ALLOC(n, ulong);
        scratch = FLINT_ARRAY_ALLOC(2 * n, ulong);

        flint_mpn_rrandom(m, state, n);
        m[0] |= 1;
        m[n - 1] |= (UWORD(1) << (FLINT_BITS - 1 -
                        n_randint(state, FLINT_BITS)));

        /* m * m^(-1) ≡ 1 mod B^n */
        _flint_mpn_binvert(v, m, n);
        flint_mpn_mullow_n(scratch, m, v, n);
        if (scratch[0] != 1 || (n > 1 && !flint_mpn_zero_p(scratch + 1, n - 1)))
            TEST_FUNCTION_FAIL("binvert: n = %wd\n", n);

        flint_mpn_copyi(minv, v, n);
        mpn_neg(minv, minv, n);

        /* X_hi < m keeps X < m * B^n */
        for (j = 0; j < 3; j++)
        {
            flint_mpn_rrandom(X, state, 2 * n);
            if (n_randint(state, 8) == 0)
                flint_mpn_zero(X, n);          /* the q = 0 branch */
            mpn_tdiv_qr(scratch, X + n, 0, X + n, n, m, n);

            /* crafted low halves drive the guard limb through its
               decision boundary: X_lo = s makes lo(q*m) = B^n - s
               (true guard UWORD_MAX-ish), X_lo = B^n - s makes
               lo(q*m) = s (true guard 0 with a wrapped approximation) */
            if (j == 1)
            {
                ulong s = 1 + n_randint(state, 2 * n + 4);
                flint_mpn_zero(X, n);
                X[0] = s;
            }
            else if (j == 2)
            {
                ulong s = 1 + n_randint(state, 2 * n + 4);
                slong i;
                for (i = 0; i < n; i++)
                    X[i] = UWORD_MAX;
                mpn_sub_1(X, X, n, s - 1);     /* X_lo = B^n - s */
            }

            _flint_mpn_redc_n(t, X, m, minv, n, scratch);

            if (!redc_ok(t, X, m, n))
                TEST_FUNCTION_FAIL("redc: n = %wd, j = %wd\n", n, j);
        }

        flint_free(m);
        flint_free(minv);
        flint_free(v);
        flint_free(X);
        flint_free(t);
        flint_free(scratch);
    }

    TEST_FUNCTION_END(state);
}
