/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fft_small.h"
#include "machine_vectors.h"

void test_mul(mpn_ctx_t R, ulong maxsize, ulong nreps, flint_rand_t state)
{
    ulong * a, * b, * c, *d;

    a = FLINT_ARRAY_ALLOC(maxsize, ulong);
    b = FLINT_ARRAY_ALLOC(maxsize, ulong);
    c = FLINT_ARRAY_ALLOC(maxsize, ulong);
    d = FLINT_ARRAY_ALLOC(maxsize, ulong);

    for (ulong rep = 0; rep < nreps; rep++)
    {
        ulong an = 2 + n_randint(state, maxsize - 4);
        ulong bn = 1 + n_randint(state, n_min(an, maxsize - an));

        for (ulong i = 0; i < maxsize; i++)
        {
            a[i] = n_randlimb(state);
            b[i] = n_randlimb(state);
            c[i] = n_randlimb(state);
            d[i] = n_randlimb(state);
        }

        mpn_ctx_mpn_mul(R, d, a, an, b, bn);
        mpn_mul(c, a, an, b, bn);
        for (ulong i = 0; i < an + bn; i++)
        {
            if (c[i] != d[i])
            {
                flint_printf("\nFAILED\n");
                flint_printf("an = %wu, bn = %wu\n", an, bn);
                flint_printf("limb[%wu] = 0x%wx should be 0x%wx\n", i, d[i], c[i]);
                fflush(stdout);
                flint_abort();
            }
        }

        mpn_ctx_mpn_mul(R, d, b, bn, b, bn);
        mpn_sqr(c, b, bn);
        for (ulong i = 0; i < 2 * bn; i++)
        {
            if (c[i] != d[i])
            {
                flint_printf("\nFAILED (squaring)\n");
                flint_printf("bn = %wu\n", bn);
                flint_printf("limb[%wu] = 0x%wx should be 0x%wx\n", i, d[i], c[i]);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_set_num_threads(1 + n_randint(state, 10));
    }

    flint_free(a);
    flint_free(b);
    flint_free(c);
    flint_free(d);
}

TEST_FUNCTION_START(mpn_ctx_mpn_mul, state)
{
    {
        mpn_ctx_t R;
        mpn_ctx_init(R, UWORD(0x0003f00000000001));
        test_mul(R, 10000, 1000 * flint_test_multiplier(), state);
        test_mul(R, 50000, 100 * flint_test_multiplier(), state);
        mpn_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}
