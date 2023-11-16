/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sin_cos_taylor_rs, state)
{
    slong iter;

    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        mp_ptr x, y1s, y1c, y2s, y2c, t;
        mp_limb_t err1, err2;
        ulong N;
        mp_size_t xn;
        int cmp, result;

        N = n_randint(state, 144 - 1);
        xn = 1 + n_randint(state, 20);

        x = flint_malloc(sizeof(mp_limb_t) * xn);
        y1s = flint_malloc(sizeof(mp_limb_t) * xn);
        y1c = flint_malloc(sizeof(mp_limb_t) * xn);
        y2s = flint_malloc(sizeof(mp_limb_t) * xn);
        y2c = flint_malloc(sizeof(mp_limb_t) * xn);
        t = flint_malloc(sizeof(mp_limb_t) * xn);

        flint_mpn_rrandom(x, state->gmp_state, xn);
        flint_mpn_rrandom(y1s, state->gmp_state, xn);
        flint_mpn_rrandom(y1c, state->gmp_state, xn);
        flint_mpn_rrandom(y2s, state->gmp_state, xn);
        flint_mpn_rrandom(y2c, state->gmp_state, xn);
        x[xn - 1] &= (LIMB_ONES >> 4);

        _arb_sin_cos_taylor_naive(y1s, y1c, &err1, x, xn, N);
        _arb_sin_cos_taylor_rs(y2s, y2c, &err2, x, xn, N, 0, 1);

        cmp = mpn_cmp(y1s, y2s, xn);

        if (cmp == 0)
        {
            result = 1;
        }
        else if (cmp > 0)
        {
            mpn_sub_n(t, y1s, y2s, xn);
            result = flint_mpn_zero_p(t + 1, xn - 1) && (t[0] <= err2);
        }
        else
        {
            mpn_sub_n(t, y2s, y1s, xn);
            result = flint_mpn_zero_p(t + 1, xn - 1) && (t[0] <= err2);
        }

        if (!result)
        {
            flint_printf("FAIL\n");
            flint_printf("N = %wd xn = %wd\n", N, xn);
            flint_printf("x =");
            flint_mpn_debug(x, xn);
            flint_printf("y1s =");
            flint_mpn_debug(y1s, xn);
            flint_printf("y2s =");
            flint_mpn_debug(y2s, xn);
            flint_abort();
        }

        cmp = mpn_cmp(y1c, y2c, xn);

        if (cmp == 0)
        {
            result = 1;
        }
        else if (cmp > 0)
        {
            mpn_sub_n(t, y1c, y2c, xn);
            result = flint_mpn_zero_p(t + 1, xn - 1) && (t[0] <= err2);
        }
        else
        {
            mpn_sub_n(t, y2c, y1c, xn);
            result = flint_mpn_zero_p(t + 1, xn - 1) && (t[0] <= err2);
        }

        if (!result)
        {
            flint_printf("FAIL\n");
            flint_printf("N = %wd xn = %wd\n", N, xn);
            flint_printf("x =");
            flint_mpn_debug(x, xn);
            flint_printf("y1c =");
            flint_mpn_debug(y1c, xn);
            flint_printf("y2c =");
            flint_mpn_debug(y2c, xn);
            flint_abort();
        }

        flint_free(x);
        flint_free(y1s);
        flint_free(y1c);
        flint_free(y2s);
        flint_free(y2c);
        flint_free(t);
    }

    TEST_FUNCTION_END(state);
}
