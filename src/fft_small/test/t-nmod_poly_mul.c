/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "fft_small.h"

TEST_FUNCTION_START(_nmod_poly_mul_mid_mpn_ctx, state)
{
    flint_bitcnt_t nbits;
    mpn_ctx_t R;
    nmod_t mod;

    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    /* (slow) test bug where 3 instead of 4 primes were used */
#if 0
    {
        ulong * a, *b, * c, * d;
        ulong an, zn, zl, zh, sz, i;

        nbits = FLINT_BITS;

        {
            nmod_init(&mod, UWORD(18446744073709551557));

            an = 25000000;
            zn = an + an - 1;
            zl = n_randint(state, zn+10);
            zh = n_randint(state, zn+20);

            sz = FLINT_MAX(zl, zh);
            sz = FLINT_MAX(sz, zn);

            a = FLINT_ARRAY_ALLOC(an, ulong);
            b = FLINT_ARRAY_ALLOC(an, ulong);
            c = FLINT_ARRAY_ALLOC(sz, ulong);
            d = FLINT_ARRAY_ALLOC(sz, ulong);

            for (i = 0; i < an; i++)
                b[i] = a[i] = n_randint(state, mod.n);

            flint_mpn_zero(c, sz);
            _nmod_poly_mul_KS(c, a, an, a, an, 0, mod);
            _nmod_poly_mul_mid_mpn_ctx(d, zl, zh, a, an, b, an, mod, R);

            for (i = zl; i < zh; i++)
            {
                if (c[i] != d[i-zl])
                {
                    flint_printf("(huge) mulmid error at index %wu\n", i);
                    flint_printf("zl=%wu, zh=%wu, an=%wu\n", zl, zh, an);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            flint_free(a);
            flint_free(b);
            flint_free(c);
            flint_free(d);
        }
    }
#endif

    /* Check squaring */
    {
        ulong * a, *b, * c, * d;
        ulong an, zn, zl, zh, sz, i, reps;

        for (reps = 0; reps < 1000 * flint_test_multiplier(); reps++)
        {
            flint_set_num_threads(1 + n_randint(state, 10));

            /* 1 <= nbits <= FLINT_BITS */
            nbits = 1 + n_randint(state, FLINT_BITS);
            nmod_init(&mod, n_randbits(state, nbits));

            an = 1 + n_randint(state, 7000);
            zn = an + an - 1;
            zl = n_randint(state, zn+10);
            zh = n_randint(state, zn+20);

            sz = FLINT_MAX(zl, zh);
            sz = FLINT_MAX(sz, zn);

            a = FLINT_ARRAY_ALLOC(an, ulong);
            b = FLINT_ARRAY_ALLOC(an, ulong);
            c = FLINT_ARRAY_ALLOC(sz, ulong);
            d = FLINT_ARRAY_ALLOC(sz, ulong);

            for (i = 0; i < an; i++)
                b[i] = a[i] = n_randint(state, mod.n);

            flint_mpn_zero(c, sz);
            _nmod_poly_mul_KS(c, a, an, a, an, 0, mod);
            _nmod_poly_mul_mid_mpn_ctx(d, zl, zh, a, an, b, an, mod, R);

            for (i = zl; i < zh; i++)
            {
                if (c[i] != d[i-zl])
                {
                    flint_printf("(squaring) mulmid error at index %wu\n", i);
                    flint_printf("zl=%wu, zh=%wu, an=%wu\n", zl, zh, an);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            flint_free(a);
            flint_free(b);
            flint_free(c);
            flint_free(d);
        }
    }

    /* Check multiplication */
    {
        ulong * a, * b, * c, * d;
        ulong an, bn, zn, zl, zh, sz, i, reps;

        for (reps = 0; reps < 1000 * flint_test_multiplier(); reps++)
        {
            flint_set_num_threads(1 + n_randint(state, 10));

            /* 1 <= nbits <= FLINT_BITS */
            nbits = 1 + n_randint(state, FLINT_BITS);
            nmod_init(&mod, n_randbits(state, nbits));

            an = 1 + n_randint(state, 7000);
            bn = 1 + n_randint(state, an);
            zn = an + bn - 1;
            zl = n_randint(state, zn+10);
            zh = n_randint(state, zn+20);

            sz = FLINT_MAX(zl, zh);
            sz = FLINT_MAX(sz, zn);

            a = FLINT_ARRAY_ALLOC(an, ulong);
            b = FLINT_ARRAY_ALLOC(bn, ulong);
            c = FLINT_ARRAY_ALLOC(sz, ulong);
            d = FLINT_ARRAY_ALLOC(sz, ulong);

            for (i = 0; i < an; i++)
                a[i] = n_randint(state, mod.n);

            for (i = 0; i < bn; i++)
                b[i] = n_randint(state, mod.n);

            flint_mpn_zero(c, sz);
            _nmod_poly_mul_KS(c, a, an, b, bn, 0, mod);
            _nmod_poly_mul_mid_mpn_ctx(d, zl, zh, a, an, b, bn, mod, R);

            for (i = zl; i < zh; i++)
            {
                if (c[i] != d[i-zl])
                {
                    flint_printf("mulmid error at index %wu\n", i);
                    flint_printf("zl=%wu, zh=%wu, an=%wu, bn=%wu\n", zl, zh, an, bn);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            flint_free(a);
            flint_free(b);
            flint_free(c);
            flint_free(d);
        }
    }

    mpn_ctx_clear(R);

    TEST_FUNCTION_END(state);
}
