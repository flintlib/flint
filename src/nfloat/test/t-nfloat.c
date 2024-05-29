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
#include "acf.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "nfloat.h"

TEST_FUNCTION_START(nfloat, state)
{
    gr_ctx_t ctx;
    gr_ctx_t ctx2;
    slong prec;
    slong iter;

    for (prec = NFLOAT_MIN_LIMBS * FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        nfloat_ctx_init(ctx, prec, 0);

        gr_test_floating_point(ctx, 100 * flint_test_multiplier(), 0);

        {
            gr_ptr x, x2;
            slong i;

            gr_ctx_init_real_arb(ctx2, 2 + n_randint(state, 500));
            x = gr_heap_init(ctx);
            x2 = gr_heap_init(ctx2);
            for (i = 0; i < 10; i++)
            {
                GR_MUST_SUCCEED(gr_randtest(x, state, ctx));
                GR_MUST_SUCCEED(gr_set_other(x2, x, ctx, ctx2));
                GR_MUST_SUCCEED(gr_set_other(x, x2, ctx2, ctx));
            }
            gr_heap_clear(x, ctx);
            gr_heap_clear(x2, ctx2);
            gr_ctx_clear(ctx2);

            gr_ctx_init_fmpz(ctx2);
            x = gr_heap_init(ctx);
            x2 = gr_heap_init(ctx2);
            for (i = 0; i < 10; i++)
            {
                GR_MUST_SUCCEED(gr_randtest(x2, state, ctx2));
                GR_MUST_SUCCEED(gr_set_other(x, x2, ctx2, ctx));
            }
            gr_heap_clear(x, ctx);
            gr_heap_clear(x2, ctx2);
            gr_ctx_clear(ctx2);

            gr_ctx_init_fmpq(ctx2);
            x = gr_heap_init(ctx);
            x2 = gr_heap_init(ctx2);
            for (i = 0; i < 10; i++)
            {
                GR_MUST_SUCCEED(gr_randtest(x2, state, ctx2));
                GR_MUST_SUCCEED(gr_set_other(x, x2, ctx2, ctx));
            }
            gr_heap_clear(x, ctx);
            gr_heap_clear(x2, ctx2);
            gr_ctx_clear(ctx2);

            gr_ctx_init_real_qqbar(ctx2);
            x = gr_heap_init(ctx);
            x2 = gr_heap_init(ctx2);
            for (i = 0; i < 3; i++)
            {
                GR_MUST_SUCCEED(gr_randtest(x2, state, ctx2));
                GR_MUST_SUCCEED(gr_set_other(x, x2, ctx2, ctx));
            }
            gr_heap_clear(x, ctx);
            gr_heap_clear(x2, ctx2);
            gr_ctx_clear(ctx2);
        }

        {
            gr_ptr tol1, tol;
            slong i, reps;

            gr_ctx_init_real_arb(ctx2, prec + 64);

            tol1 = gr_heap_init(ctx);
            tol = gr_heap_init(ctx2);

            GR_IGNORE(gr_one(tol, ctx2));
            GR_IGNORE(gr_mul_2exp_si(tol, tol, -prec + 2, ctx2));

            reps = (prec <= 256 ? 10000 : 1) * flint_test_multiplier();

            for (i = 0; i < reps; i++)
            {
                gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_add, ctx2, tol, state, 0);
                gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_sub, ctx2, tol, state, 0);
            }

            reps = (prec <= 256 ? 10000 : 100) * flint_test_multiplier();
            for (i = 0; i < reps; i++)
                gr_test_approx_dot(ctx, ctx2, 30, tol, state, 0);

            reps = (prec <= 256 ? 1000 : 1) * flint_test_multiplier();

            for (i = 0; i < reps; i++)
            {
                gr_test_cmp_fun(ctx, (gr_method_binary_op_get_int) gr_cmp, ctx2, state, 0);
                gr_test_cmp_fun(ctx, (gr_method_binary_op_get_int) gr_cmpabs, ctx2, state, 0);

                gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_mul, ctx2, tol, state, 0);
                gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_div, ctx2, tol, state, 0);
                gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_pow, ctx2, tol, state, 0);

                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_neg, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_abs, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_sgn, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_sqr, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_inv, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_sqrt, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_rsqrt, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_floor, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_ceil, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_nint, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_trunc, ctx2, tol, state, 0);

                GR_IGNORE(gr_one(tol1, ctx));
                GR_IGNORE(gr_mul_2exp_si(tol1, tol1, -prec + 4, ctx));
                gr_test_approx_binary_op_type_variants(ctx,
                    "mul",
                    (gr_method_binary_op) gr_mul,
                    (gr_method_binary_op_ui) gr_mul_ui,
                    (gr_method_binary_op_si) gr_mul_si,
                    (gr_method_binary_op_fmpz) gr_mul_fmpz,
                    (gr_method_binary_op_fmpq) gr_mul_fmpq,
                    0, 0, tol1, state, 0);

                gr_test_approx_binary_op_type_variants(ctx,
                    "div",
                    (gr_method_binary_op) gr_div,
                    (gr_method_binary_op_ui) gr_div_ui,
                    (gr_method_binary_op_si) gr_div_si,
                    (gr_method_binary_op_fmpz) gr_div_fmpz,
                    (gr_method_binary_op_fmpq) gr_div_fmpq,
                    0, 0, tol1, state, 0);

                /* todo: test with large values also */
                GR_IGNORE(gr_one(tol1, ctx));
                GR_IGNORE(gr_mul_2exp_si(tol1, tol1, -prec + 6, ctx));
                gr_test_approx_binary_op_type_variants(ctx,
                    "pow",
                    (gr_method_binary_op) gr_pow,
                    (gr_method_binary_op_ui) gr_pow_ui,
                    (gr_method_binary_op_si) gr_pow_si,
                    (gr_method_binary_op_fmpz) gr_pow_fmpz,
                    (gr_method_binary_op_fmpq) gr_pow_fmpq,
                    0, 1, tol1, state, 0);
            }

            reps = (prec <= 256 ? 10 : 0) * flint_test_multiplier();

            for (i = 0; i < reps; i++)
            {
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_exp, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_expm1, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_log, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_log1p, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_sin, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_cos, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_tan, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_sinh, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_cosh, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_tanh, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_atan, ctx2, tol, state, 0);
            }

            reps = (prec <= 256 ? 1 : 0) * flint_test_multiplier();

            for (i = 0; i < reps; i++)
            {
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_gamma, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_zeta, ctx2, tol, state, 0);
            }

            gr_heap_clear(tol1, ctx);
            gr_heap_clear(tol, ctx2);
            gr_ctx_clear(ctx2);
        }

        gr_ctx_clear(ctx);
    }

    nfloat_ctx_init(ctx, FLINT_BITS, 0);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        ulong x[3];
        ulong r[NFLOAT_MAX_ALLOC];
        ulong s[NFLOAT_MAX_ALLOC];
        int s1, s2;
        slong exp;
        int sgn;

        x[0] = n_randint(state, 2) ? 0 : n_randtest(state);
        x[1] = n_randint(state, 2) ? 0 : n_randtest(state);
        x[2] = n_randint(state, 2) ? 0 : n_randtest(state);
        exp = (slong) n_randint(state, 100) - 100;
        sgn = n_randint(state, 2);

        s1 = nfloat_1_set_3_2exp(r, x[2], x[1], x[0], exp, sgn, ctx);
        s2 = nfloat_set_mpn_2exp(s, x, 3, exp, sgn, ctx);

        if (s1 != s2 || nfloat_equal(r, s, ctx) == T_FALSE)
        {
            flint_printf("FAIL: nfloat_1_set_3_2exp\n");
            flint_mpn_debug(x, 3);
            gr_println(r, ctx);
            gr_println(s, ctx);
            flint_abort();
        }
    }

    gr_ctx_clear(ctx);

    nfloat_ctx_init(ctx, 2 * FLINT_BITS, 0);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        ulong x[4];
        ulong r[NFLOAT_MAX_ALLOC];
        ulong s[NFLOAT_MAX_ALLOC];
        int s1, s2;
        slong exp;
        int sgn;

        x[0] = n_randint(state, 2) ? 0 : n_randtest(state);
        x[1] = n_randint(state, 2) ? 0 : n_randtest(state);
        x[2] = n_randint(state, 2) ? 0 : n_randtest(state);
        x[3] = n_randint(state, 2) ? 0 : n_randtest(state);
        exp = (slong) n_randint(state, 100) - 100;
        sgn = n_randint(state, 2);

        s1 = nfloat_2_set_4_2exp(r, x[3], x[2], x[1], x[0], exp, sgn, ctx);
        s2 = nfloat_set_mpn_2exp(s, x, 4, exp, sgn, ctx);

        if (s1 != s2 || nfloat_equal(r, s, ctx) == T_FALSE)
        {
            flint_printf("FAIL: nfloat_2_set_4_2exp\n");
            flint_mpn_debug(x, 4);
            gr_println(r, ctx);
            gr_println(s, ctx);
            flint_abort();
        }
    }

    gr_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
