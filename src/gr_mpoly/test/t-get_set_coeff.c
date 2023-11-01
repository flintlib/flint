/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_get_set_coeff, state)
{
    slong i, j;

    for (i = 0; i < 100; i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f;
        gr_ptr c, d;
        slong len, k;
        flint_bitcnt_t exp_bits;
        int status;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 20);

        gr_mpoly_init(f, mctx, cctx);
        c = gr_heap_init(cctx);
        d = gr_heap_init(cctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;

        GR_MUST_SUCCEED(gr_mpoly_randtest_bits(f, state, len, exp_bits, mctx, cctx));

        for (j = 0; j < 10; j++)
        {
            ulong * exp = FLINT_ARRAY_ALLOC(mctx->nvars, ulong);

            status = GR_SUCCESS;

            GR_MUST_SUCCEED(gr_randtest(c, state, cctx));
            for (k = 0; k < mctx->nvars; k++)
                exp[k] = n_randtest(state);

            status |= gr_mpoly_set_coeff_scalar_ui(f, c, exp, mctx, cctx);
            gr_mpoly_assert_canonical(f, mctx, cctx);
            status |= gr_mpoly_get_coeff_scalar_ui(d, f, exp, mctx, cctx);

            if (status == GR_SUCCESS && gr_equal(c, d, cctx) == T_FALSE)
            {
                flint_printf("FAIL: scalar_ui\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("c = "); gr_print(c, cctx); flint_printf("\n");
                flint_printf("d = "); gr_print(d, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            flint_free(exp);
        }

        gr_mpoly_clear(f, mctx, cctx);
        gr_heap_clear(c, cctx);
        gr_heap_clear(d, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

    for (i = 0; i < 100; i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f;
        gr_ptr c, d;
        slong len;
        flint_bitcnt_t exp_bits;
        int status;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 20);

        gr_mpoly_init(f, mctx, cctx);
        c = gr_heap_init(cctx);
        d = gr_heap_init(cctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;

        GR_MUST_SUCCEED(gr_mpoly_randtest_bits(f, state, len, exp_bits, mctx, cctx));

        for (j = 0; j < 10; j++)
        {
            fmpz * exp = _fmpz_vec_init(mctx->nvars);

            status = GR_SUCCESS;

            GR_MUST_SUCCEED(gr_randtest(c, state, cctx));
            _fmpz_vec_randtest(exp, state, mctx->nvars, exp_bits);

            status |= gr_mpoly_set_coeff_scalar_fmpz(f, c, exp, mctx, cctx);
            gr_mpoly_assert_canonical(f, mctx, cctx);
            status |= gr_mpoly_get_coeff_scalar_fmpz(d, f, exp, mctx, cctx);

            if (status == GR_SUCCESS && gr_equal(c, d, cctx) == T_FALSE)
            {
                flint_printf("FAIL: scalar_fmpz\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("c = "); gr_print(c, cctx); flint_printf("\n");
                flint_printf("d = "); gr_print(d, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            _fmpz_vec_clear(exp, mctx->nvars);
        }

        gr_mpoly_clear(f, mctx, cctx);
        gr_heap_clear(c, cctx);
        gr_heap_clear(d, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
