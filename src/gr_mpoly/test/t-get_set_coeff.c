/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_get_set_coeff, state)
{
    slong i, j;

    for (i = 0; i < 100; i++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f;
        gr_ptr c, d;
        slong len, k;
        flint_bitcnt_t exp_bits;
        int status;
        slong nvars;
        int which;

        gr_ctx_init_random(cctx, state);
        gr_mpoly_ctx_init_rand(ctx, state, cctx, 20);
        nvars = GR_MPOLY_NVARS(ctx);

        gr_mpoly_init(f, ctx);
        c = gr_heap_init(cctx);
        d = gr_heap_init(cctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;

        GR_MUST_SUCCEED(gr_mpoly_randtest_bits(f, state, len, exp_bits, ctx));

        for (j = 0; j < 10; j++)
        {
            ulong * exp = FLINT_ARRAY_ALLOC(nvars, ulong);

            status = GR_SUCCESS;

            for (k = 0; k < nvars; k++)
                exp[k] = n_randtest(state);

            which = n_randint(state, 5);

            if (which == 0)
            {
                GR_MUST_SUCCEED(gr_randtest(c, state, cctx));
                status |= gr_mpoly_set_coeff_scalar_ui(f, c, exp, ctx);
            }
            else if (which == 1)
            {
                ulong uc = n_randtest(state);
                status |= gr_set_ui(c, uc, cctx);
                status |= gr_mpoly_set_coeff_ui_ui(f, uc, exp, ctx);
            }
            else if (which == 2)
            {
                slong sc = n_randtest(state);
                status |= gr_set_si(c, sc, cctx);
                status |= gr_mpoly_set_coeff_si_ui(f, sc, exp, ctx);
            }
            else if (which == 3)
            {
                fmpz_t zc;
                fmpz_init(zc);
                fmpz_randtest(zc, state, 100);
                status |= gr_set_fmpz(c, zc, cctx);
                status |= gr_mpoly_set_coeff_fmpz_ui(f, zc, exp, ctx);
                fmpz_clear(zc);
            }
            else if (which == 4)
            {
                fmpq_t qc;
                fmpq_init(qc);
                fmpq_randtest(qc, state, 100);
                status |= gr_set_fmpq(c, qc, cctx);
                status |= gr_mpoly_set_coeff_fmpq_ui(f, qc, exp, ctx);
                fmpq_clear(qc);
            }

            gr_mpoly_assert_canonical(f, ctx);
            status |= gr_mpoly_get_coeff_scalar_ui(d, f, exp, ctx);

            if (status == GR_SUCCESS && gr_equal(c, d, cctx) == T_FALSE)
            {
                flint_printf("FAIL: scalar_ui\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("c = "); gr_print(c, cctx); flint_printf("\n");
                flint_printf("d = "); gr_print(d, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            flint_free(exp);
        }

        GR_MUST_SUCCEED(gr_mpoly_randtest_bits(f, state, len, exp_bits, ctx));

        for (j = 0; j < 10; j++)
        {
            fmpz * exp = _fmpz_vec_init(nvars);

            status = GR_SUCCESS;

            _fmpz_vec_randtest(exp, state, nvars, exp_bits);

            which = n_randint(state, 5);

            if (which == 0)
            {
                GR_MUST_SUCCEED(gr_randtest(c, state, cctx));
                status |= gr_mpoly_set_coeff_scalar_fmpz(f, c, exp, ctx);
            }
            else if (which == 1)
            {
                ulong uc = n_randtest(state);
                status |= gr_set_ui(c, uc, cctx);
                status |= gr_mpoly_set_coeff_ui_fmpz(f, uc, exp, ctx);
            }
            else if (which == 2)
            {
                slong sc = n_randtest(state);
                status |= gr_set_si(c, sc, cctx);
                status |= gr_mpoly_set_coeff_si_fmpz(f, sc, exp, ctx);
            }
            else if (which == 3)
            {
                fmpz_t zc;
                fmpz_init(zc);
                fmpz_randtest(zc, state, 100);
                status |= gr_set_fmpz(c, zc, cctx);
                status |= gr_mpoly_set_coeff_fmpz_fmpz(f, zc, exp, ctx);
                fmpz_clear(zc);
            }
            else if (which == 4)
            {
                fmpq_t qc;
                fmpq_init(qc);
                fmpq_randtest(qc, state, 100);
                status |= gr_set_fmpq(c, qc, cctx);
                status |= gr_mpoly_set_coeff_fmpq_fmpz(f, qc, exp, ctx);
                fmpq_clear(qc);
            }

            gr_mpoly_assert_canonical(f, ctx);
            status |= gr_mpoly_get_coeff_scalar_fmpz(d, f, exp, ctx);

            if (status == GR_SUCCESS && gr_equal(c, d, cctx) == T_FALSE)
            {
                flint_printf("FAIL: scalar_fmpz\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("c = "); gr_print(c, cctx); flint_printf("\n");
                flint_printf("d = "); gr_print(d, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            _fmpz_vec_clear(exp, nvars);
        }

        gr_mpoly_clear(f, ctx);
        gr_heap_clear(c, cctx);
        gr_heap_clear(d, cctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
