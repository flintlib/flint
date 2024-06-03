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
#include "gr_vec.h"
#include "gr_special.h"
#include "nfloat.h"

TEST_FUNCTION_START(nfloat_complex, state)
{
    gr_ctx_t ctx;
    gr_ctx_t ctx2;
    slong prec;
    slong iter, reps;
    gr_ptr tol1, tol;

    for (prec = NFLOAT_MIN_LIMBS * FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        nfloat_complex_ctx_init(ctx, prec, 0);
        gr_ctx_init_complex_acb(ctx2, prec + 64);

        gr_test_floating_point(ctx, 100 * flint_test_multiplier(), 0);

        tol1 = gr_heap_init(ctx);
        tol = gr_heap_init(ctx2);

        GR_IGNORE(gr_one(tol, ctx2));
        GR_IGNORE(gr_mul_2exp_si(tol, tol, -prec + 3, ctx2));
        reps = (prec <= 256 ? 10000 : 1) * flint_test_multiplier();
        for (iter = 0; iter < reps; iter++)
        {
            gr_test_cmp_fun(ctx, (gr_method_binary_op_get_int) gr_cmp, ctx2, state, 0);
            gr_test_cmp_fun(ctx, (gr_method_binary_op_get_int) gr_cmpabs, ctx2, state, 0);
            gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_add, ctx2, tol, state, 0);
            gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_sub, ctx2, tol, state, 0);
            gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_mul, ctx2, tol, state, 0);
            gr_test_approx_binary_op(ctx, (gr_method_binary_op) gr_div, ctx2, tol, state, 0);

            gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_neg, ctx2, tol, state, 0);
            gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_sqr, ctx2, tol, state, 0);
            gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_inv, ctx2, tol, state, 0);
            gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_abs, ctx2, tol, state, 0);
        }

        GR_IGNORE(gr_one(tol, ctx2));
        GR_IGNORE(gr_mul_2exp_si(tol, tol, -prec + 4, ctx2));
        reps = (prec <= 256 ? 10000 : 100) * flint_test_multiplier();
        for (iter = 0; iter < reps; iter++)
            gr_test_approx_dot(ctx, ctx2, 10, tol, state, 0);

        gr_heap_clear(tol1, ctx);
        gr_heap_clear(tol, ctx2);
        gr_ctx_clear(ctx2);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}