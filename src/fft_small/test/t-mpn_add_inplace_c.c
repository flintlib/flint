/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fft_small.h"
#include "crt_helpers.h"

#define TEST_ADD(n) \
        flint_mpn_copyi(a, d, n); \
        mpn_add_n(c, a, b, n); \
        multi_add_ ## n(a, b); \
        if (mpn_cmp(a, c, n) != 0) \
        { \
            flint_printf("FAIL: add_%i\n", n); \
            flint_abort(); \
        }

#define TEST_SUB(n) \
        flint_mpn_copyi(a, d, n); \
        mpn_sub_n(c, a, b, n); \
        multi_sub_ ## n(a, b); \
        if (mpn_cmp(a, c, n) != 0) \
        { \
            flint_printf("FAIL: sub_%i\n", n); \
            flint_abort(); \
        }

TEST_FUNCTION_START(flint_mpn_add_inplace_c, state)
{
    slong iter;

    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        mp_limb_t a[10], b[10], c[10];
        mp_size_t an, bn;
        unsigned char cf, c1, c2;

        bn = 1 + n_randint(state, 4);
        an = bn + n_randint(state, 4);

        flint_mpn_rrandom(a, state->gmp_state, an);
        flint_mpn_rrandom(b, state->gmp_state, bn);
        flint_mpn_copyi(c, a, an);
        cf = n_randint(state, 2);

        c1 = flint_mpn_add_inplace_c(a, an, b, bn, cf);

        c2 = mpn_add(c, c, an, b, bn);
        c2 += mpn_add_1(c, c, an, cf);

        if (c1 != c2 || mpn_cmp(a, c, an) != 0)
        {
            flint_printf("FAIL\n");
            flint_abort();
        }
    }

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        mp_limb_t a[8], b[8], c[8], d[8];

        flint_mpn_rrandom(a, state->gmp_state, 8);
        flint_mpn_rrandom(b, state->gmp_state, 8);
        flint_mpn_copyi(c, a, 8);
        flint_mpn_copyi(d, a, 8);

        TEST_ADD(1)
        TEST_ADD(2)
        TEST_ADD(3)
        TEST_ADD(4)
        TEST_ADD(5)
        TEST_ADD(6)
        TEST_ADD(7)
        TEST_ADD(8)

        TEST_SUB(1)
        TEST_SUB(2)
        TEST_SUB(3)
        TEST_SUB(4)
        TEST_SUB(5)
        TEST_SUB(6)
        TEST_SUB(7)
        TEST_SUB(8)
    }

    TEST_FUNCTION_END(state);
}
