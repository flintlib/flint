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

TEST_FUNCTION_START(flint_mpn_powm, state)
{
    for (slong iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        int large = n_randint(state, 12) == 0;
        int huge = n_randint(state, 40) == 0;
        mp_size_t mn = huge ? 1550 + n_randint(state, 500)
                            : 1 + n_randint(state, large ? 350 : 25);
        mp_size_t en = (large || huge) ? 1 : n_randint(state, 5);
        nn_ptr m, b, e, r, r2, q;
        mpz_t mz, bz, ez, rz, cz;
        int mode;

        m = FLINT_ARRAY_ALLOC(mn, ulong);
        b = FLINT_ARRAY_ALLOC(mn, ulong);
        e = FLINT_ARRAY_ALLOC(FLINT_MAX(en, 1), ulong);
        r = FLINT_ARRAY_ALLOC(mn, ulong);
        r2 = FLINT_ARRAY_ALLOC(mn, ulong);
        q = FLINT_ARRAY_ALLOC(mn + 1, ulong);

        /* modulus: random, sometimes even with a big 2-adic part,
           sometimes a pure power of two */
        flint_mpn_rrandom(m, state, mn);
        m[mn - 1] |= (UWORD(1) << (FLINT_BITS - 1 -
                          n_randint(state, FLINT_BITS)));
        mode = n_randint(state, 4);
        if (mode == 1)
            m[0] |= 1;                                   /* odd */
        else if (mode == 2 && mn > 1)
            flint_mpn_zero(m, n_randint(state, mn - 1) + 1),  /* even, big val2 */
            m[mn - 1] |= 1 + n_randint(state, 255);
        else if (mode == 3)
        {
            flint_mpn_zero(m, mn);                       /* m = 2^k */
            m[mn - 1] = UWORD(1) << n_randint(state, FLINT_BITS);
        }
        if (flint_mpn_zero_p(m, mn))
            m[0] = 1;
        while (m[mn - 1] == 0)
            mn--;

        /* base: random < m, or a small prime to hit the scalar path */
        if (n_randint(state, 3) == 0)
        {
            ulong small[] = { 0, 1, 2, 3, 5, 7, 10, 255 };
            flint_mpn_zero(b, mn);
            b[0] = small[n_randint(state, 8)];
            if (mn == 1 || (mn == 1 && b[0] >= m[0]))
                b[0] = b[0] % m[0];
            if (mn == 1)
                b[0] %= m[0];
        }
        else
        {
            flint_mpn_rrandom(b, state, mn);
            mpn_tdiv_qr(q, b, 0, b, mn, m, mn);
        }
        if (mpn_cmp(b, m, mn) >= 0)
            mpn_sub_n(b, b, m, mn);

        /* exponent: random, all-ones, a power of two (long runs), or
           an explicit small value (the preinvn fast path) */
        if (huge && en > 0)
        {
            e[0] = 17 + n_randint(state, 100);
        }
        else if (en > 0 && n_randint(state, 4) == 0)
        {
            ulong smalle[] = { 2, 3, 4, 5, 16, 100, 65537 };
            en = 1;
            e[0] = smalle[n_randint(state, 7)];
        }
        else if (en > 0)
        {
            slong i;
            switch (n_randint(state, 3))
            {
                case 0: flint_mpn_rrandom(e, state, en); break;
                case 1: for (i = 0; i < en; i++) e[i] = UWORD_MAX; break;
                default:
                    flint_mpn_zero(e, en);
                    e[en - 1] = UWORD(1) << n_randint(state, FLINT_BITS);
            }
        }

        flint_mpn_powm(r, b, e, en, m, mn);

        /* the preinvn entry must agree everywhere */
        {
            flint_bitcnt_t norm = flint_clz(m[mn - 1]);
            nn_ptr msh = FLINT_ARRAY_ALLOC(2 * mn, ulong);
            nn_ptr dinv = msh + mn;
            if (norm)
                mpn_lshift(msh, m, mn, norm);
            else
                flint_mpn_copyi(msh, m, mn);
            flint_mpn_preinvn(dinv, msh, mn);
            flint_mpn_powm_preinvn(r2, b, e, en, m, mn, dinv, norm);
            if (mpn_cmp(r, r2, mn) != 0)
                TEST_FUNCTION_FAIL("preinvn entry: mn = %wd, en = %wd\n",
                                   mn, en);
            flint_free(msh);
        }

        mpz_init(mz); mpz_init(bz); mpz_init(ez); mpz_init(rz); mpz_init(cz);
        mpz_import(mz, mn, -1, sizeof(ulong), 0, 0, m);
        mpz_import(bz, mn, -1, sizeof(ulong), 0, 0, b);
        if (en > 0)
            mpz_import(ez, en, -1, sizeof(ulong), 0, 0, e);
        mpz_powm(rz, bz, ez, mz);
        mpz_import(cz, mn, -1, sizeof(ulong), 0, 0, r);

        if (mpz_cmp(rz, cz) != 0)
            TEST_FUNCTION_FAIL("powm vs mpz: mn = %wd, en = %wd, "
                               "mode = %d, b0 = %wu\n", mn, en, mode, b[0]);

        /* the stage functions, under their preconditions */
        if (en > 0 && !flint_mpn_zero_p(b, mn) && e[en - 1] != 0 &&
            mn >= 2 && !(mn == 1 && m[0] == 1))
        {
            _flint_mpn_powm_basecase(r2, b, e, en, m, mn);
            if (mpn_cmp(r, r2, mn) != 0)
                TEST_FUNCTION_FAIL("basecase: mn = %wd, en = %wd, mode = %d\n",
                                   mn, en, mode);

            if (m[0] & 1)
            {
                _flint_mpn_powm_redc(r2, b, e, en, m, mn);
                if (mpn_cmp(r, r2, mn) != 0)
                    TEST_FUNCTION_FAIL("redc: mn = %wd, en = %wd\n", mn, en);

            }
        }

        mpz_clear(mz); mpz_clear(bz); mpz_clear(ez);
        mpz_clear(rz); mpz_clear(cz);
        flint_free(m); flint_free(b); flint_free(e);
        flint_free(r); flint_free(r2); flint_free(q);
    }

    TEST_FUNCTION_END(state);
}
