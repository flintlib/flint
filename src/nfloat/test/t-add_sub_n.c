/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arf.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "nfloat.h"

TEST_FUNCTION_START(add_sub_n, state)
{
    gr_ctx_t ctx;
    slong iter, prec;
    gr_ptr a, b, r1, r2;

    /* test that add_X and sub_X specializations implement the same
       algorithm as add_n and sub_n */
    for (prec = FLINT_BITS; prec <= 4 * FLINT_BITS; prec += FLINT_BITS)
    {
        nfloat_ctx_init(ctx, prec, 0);
        GR_TMP_INIT4(a, b, r1, r2, ctx);

        for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
        {
            GR_IGNORE(gr_randtest_not_zero(a, state, ctx));
            GR_IGNORE(gr_randtest_not_zero(b, state, ctx));

            if (NFLOAT_EXP(a) < NFLOAT_EXP(b))
                nfloat_swap(a, b, ctx);

            if (n_randint(state, 2))
                NFLOAT_EXP(b) = NFLOAT_EXP(a) - 1;

            /* trigger special cases with increased probability */
            if (n_randint(state, 100) == 0)
                flint_mpn_store(NFLOAT_D(b), prec / FLINT_BITS, ~UWORD(0));

            GR_MUST_SUCCEED(_nfloat_add_n(r1, NFLOAT_D(a), NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b), NFLOAT_EXP(a) - NFLOAT_EXP(b), prec / FLINT_BITS, ctx));

            if (prec == FLINT_BITS)
                GR_MUST_SUCCEED(_nfloat_add_1(r2, NFLOAT_D(a)[0], NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b)[0], NFLOAT_EXP(a) - NFLOAT_EXP(b), ctx));
            else if (prec == 2 * FLINT_BITS)
                GR_MUST_SUCCEED(_nfloat_add_2(r2, NFLOAT_D(a), NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b), NFLOAT_EXP(a) - NFLOAT_EXP(b), ctx));
            else if (prec == 3 * FLINT_BITS)
                GR_MUST_SUCCEED(_nfloat_add_3(r2, NFLOAT_D(a), NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b), NFLOAT_EXP(a) - NFLOAT_EXP(b), ctx));
            else if (prec == 4 * FLINT_BITS)
                GR_MUST_SUCCEED(_nfloat_add_4(r2, NFLOAT_D(a), NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b), NFLOAT_EXP(a) - NFLOAT_EXP(b), ctx));

            if (gr_equal(r1, r2, ctx) != T_TRUE)
            {
                flint_printf("FAIL: add_%wd\n", prec / FLINT_BITS);
                flint_printf("a = "); gr_println(a, ctx);
                flint_printf("b = "); gr_println(b, ctx);
                flint_printf("r1 = "); gr_println(r1, ctx);
                flint_printf("r2 = "); gr_println(r2, ctx);
                flint_abort();
            }

            GR_MUST_SUCCEED(_nfloat_sub_n(r1, NFLOAT_D(a), NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b), NFLOAT_EXP(a) - NFLOAT_EXP(b), prec / FLINT_BITS, ctx));

            if (prec == FLINT_BITS)
                GR_MUST_SUCCEED(_nfloat_sub_1(r2, NFLOAT_D(a)[0], NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b)[0], NFLOAT_EXP(a) - NFLOAT_EXP(b), ctx));
            else if (prec == 2 * FLINT_BITS)
                GR_MUST_SUCCEED(_nfloat_sub_2(r2, NFLOAT_D(a), NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b), NFLOAT_EXP(a) - NFLOAT_EXP(b), ctx));
            else if (prec == 3 * FLINT_BITS)
                GR_MUST_SUCCEED(_nfloat_sub_3(r2, NFLOAT_D(a), NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b), NFLOAT_EXP(a) - NFLOAT_EXP(b), ctx));
            else if (prec == 4 * FLINT_BITS)
                GR_MUST_SUCCEED(_nfloat_sub_4(r2, NFLOAT_D(a), NFLOAT_EXP(a), NFLOAT_SGNBIT(a), NFLOAT_D(b), NFLOAT_EXP(a) - NFLOAT_EXP(b), ctx));

            if (gr_equal(r1, r2, ctx) != T_TRUE)
            {
                flint_printf("FAIL: sub_%wd\n", prec / FLINT_BITS);
                flint_printf("a = "); gr_println(a, ctx);
                flint_printf("b = "); gr_println(b, ctx);
                flint_printf("r1 = "); gr_println(r1, ctx);
                flint_printf("r2 = "); gr_println(r2, ctx);
                flint_abort();
            }
        }

        GR_TMP_CLEAR4(a, b, r1, r2, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
