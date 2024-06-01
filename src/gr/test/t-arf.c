/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"
#include "gr_special.h"

TEST_FUNCTION_START(gr_arf, state)
{
    gr_ctx_t ctx, ctx2;
    int flags = 0;
    slong prec;

    for (prec = 64; prec <= 256; prec *= 2)
    {
        gr_ctx_init_real_float_arf(ctx, prec);
        gr_test_floating_point(ctx, 100, flags);

        {
            gr_ptr tol1, tol;
            slong i, reps;

            gr_ctx_init_real_arb(ctx2, prec + 64);

            tol1 = gr_heap_init(ctx);
            tol = gr_heap_init(ctx2);

            GR_IGNORE(gr_one(tol, ctx2));
            GR_IGNORE(gr_mul_2exp_si(tol, tol, -prec + 2, ctx2));

            reps = 100 * flint_test_multiplier();

            for (i = 0; i < reps; i++)
            {
                gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_add, ctx2, tol, state, 0);
                gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_sub, ctx2, tol, state, 0);

                gr_test_approx_dot(ctx, ctx2, 10, tol, state, 0);

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
                    "add",
                    (gr_method_binary_op) gr_add,
                    (gr_method_binary_op_ui) gr_add_ui,
                    (gr_method_binary_op_si) gr_add_si,
                    (gr_method_binary_op_fmpz) gr_add_fmpz,
                    (gr_method_binary_op_fmpq) gr_add_fmpq,
                    0, 0, tol1, state, 0);

                gr_test_approx_binary_op_type_variants(ctx,
                    "sub",
                    (gr_method_binary_op) gr_sub,
                    (gr_method_binary_op_ui) gr_sub_ui,
                    (gr_method_binary_op_si) gr_sub_si,
                    (gr_method_binary_op_fmpz) gr_sub_fmpz,
                    (gr_method_binary_op_fmpq) gr_sub_fmpq,
                    0, 0, tol1, state, 0);

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

            reps = 3 * flint_test_multiplier();

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
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_gamma, ctx2, tol, state, 0);
                gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_zeta, ctx2, tol, state, 0);
            }

            gr_heap_clear(tol1, ctx);
            gr_heap_clear(tol, ctx2);
            gr_ctx_clear(ctx2);
        }

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
