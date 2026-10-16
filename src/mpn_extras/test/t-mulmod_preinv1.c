/*
    Copyright (C) 2013 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_mulmod_preinv1, state)
{
    int i;
    mpz_t a, b, d, r1;
    gmp_randstate_t st;
    mp_limb_t dinv;
    mp_size_t size;
    flint_bitcnt_t norm;

    mpz_init(a);
    mpz_init(b);
    mpz_init(d);
    mpz_init(r1);

    gmp_randinit_default(st);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_ptr A, B, D, R, R1;
        int alias;

        size = n_randint(state, 200) + 1;

        /* d of exactly size limbs, normalised or not */
        mpz_rrandomb(d, st, size * FLINT_BITS - n_randint(state, FLINT_BITS));
        mpz_rrandomb(a, st, n_randint(state, size * FLINT_BITS + 1));
        mpz_rrandomb(b, st, n_randint(state, size * FLINT_BITS + 1));

        /* reduce a, b mod d */
        mpz_fdiv_r(a, a, d);
        mpz_fdiv_r(b, b, d);

        alias = (n_randint(state, 4) == 0);
        if (alias)
            mpz_set(b, a);

        mpz_mul(r1, a, b);
        mpz_fdiv_r(r1, r1, d);

        /* normalise; the result has the same shift */
        norm = flint_clz(d->_mp_d[d->_mp_size - 1]);
        mpz_mul_2exp(a, a, norm);
        mpz_mul_2exp(b, b, norm);
        mpz_mul_2exp(d, d, norm);
        mpz_mul_2exp(r1, r1, norm);

        /* zero-padded operands of size limbs */
        A = flint_calloc(size, sizeof(mp_limb_t));
        B = flint_calloc(size, sizeof(mp_limb_t));
        D = flint_calloc(size, sizeof(mp_limb_t));
        R = flint_calloc(size, sizeof(mp_limb_t));
        R1 = flint_calloc(size, sizeof(mp_limb_t));
        mpz_export(A, NULL, -1, sizeof(mp_limb_t), 0, 0, a);
        mpz_export(B, NULL, -1, sizeof(mp_limb_t), 0, 0, b);
        mpz_export(D, NULL, -1, sizeof(mp_limb_t), 0, 0, d);
        mpz_export(R1, NULL, -1, sizeof(mp_limb_t), 0, 0, r1);

        /* for one limb, the 3/2 inverse of (d, 0) */
        dinv = flint_mpn_preinv1(D[size - 1], (size > 1) ? D[size - 2] : 0);

        flint_mpn_mulmod_preinv1(R, A, alias ? A : B, size, D, dinv, norm);

        if (mpn_cmp(R, R1, size) != 0)
            TEST_FUNCTION_FAIL("size = %wd, norm = %wu, alias = %d\n"
                "a = %{mpz}\nb = %{mpz}\nd = %{mpz}\nr1 = %{mpz}\nr = %{ulong*}\n",
                size, norm, alias, a, b, d, r1, R, size);

        flint_free(A);
        flint_free(B);
        flint_free(D);
        flint_free(R);
        flint_free(R1);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(d);
    mpz_clear(r1);

    gmp_randclear(st);

    TEST_FUNCTION_END(state);
}
