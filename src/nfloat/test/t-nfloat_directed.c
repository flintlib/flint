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

static int _arf_sqr(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
{
    return arf_mul(res, x, x, prec, rnd);
}

static int _arf_inv(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
{
    return arf_ui_div(res, 1, x, prec, rnd);
}

TEST_FUNCTION_START(nfloat_directed, state)
{
    gr_ctx_t ctx;
    gr_ctx_t ctx2;
    slong prec;
    int direction;

    for (direction = 0; direction <= 1; direction++)
    {
        for (prec = NFLOAT_MIN_LIMBS * FLINT_BITS; prec <= 320; prec += FLINT_BITS)
        {
            nfloat_ctx_init(ctx, prec, direction ? NFLOAT_RND_FLOOR : NFLOAT_RND_CEIL);

            gr_test_floating_point(ctx, 100 * flint_test_multiplier(), 0);

            /* Todo: make directed rounding tests generic */
            /* Todo: test vector operations */
            {
                slong i;
                gr_ptr x, y, z;
                arf_t ax, ay, az, az2;
                arf_rnd_t rnd = direction ? ARF_RND_FLOOR : ARF_RND_CEIL;

                x = gr_heap_init(ctx);
                y = gr_heap_init(ctx);
                z = gr_heap_init(ctx);

                arf_init(ax);
                arf_init(ay);
                arf_init(az);
                arf_init(az2);

                for (i = 0; i < 100 * flint_test_multiplier(); i++)
                {
                    GR_MUST_SUCCEED(gr_randtest(x, state, ctx));
                    GR_MUST_SUCCEED(gr_randtest(y, state, ctx));
                    GR_MUST_SUCCEED(gr_randtest(z, state, ctx));

                    nfloat_get_arf(ax, x, ctx);
                    nfloat_get_arf(ay, y, ctx);
                    nfloat_get_arf(az, z, ctx);

#define TEST_OP(gr_op, arf_op, opname) \
                    if (gr_op(z, x, y, ctx) == GR_SUCCESS) \
                    { \
                        nfloat_get_arf(az2, z, ctx); \
                        arf_op(az, ax, ay, prec, rnd); \
                        if (arf_cmp(az2, az) == ((rnd == ARF_RND_FLOOR) ? 1 : -1)) \
                        { \
                            flint_printf("FAIL: %s\n", opname); \
                            flint_printf("prec = %wd, direction = %s\n", prec, (rnd == ARF_RND_FLOOR) ? "floor" : "ceil"); \
                            flint_printf("x = %{gr}\n", x, ctx); \
                            flint_printf("y = %{gr}\n", y, ctx); \
                            flint_printf("z = %{gr}\n", z, ctx); \
                            flint_printf("ax = "); arf_printd(ax, prec / 3.32 + 1); flint_printf("\n"); \
                            flint_printf("ay = "); arf_printd(ay, prec / 3.32 + 1); flint_printf("\n"); \
                            flint_printf("az = "); arf_printd(az, prec / 3.32 + 1); flint_printf("\n"); \
                            flint_printf("az2 = "); arf_printd(az2, prec / 3.32 + 1); flint_printf("\n"); \
                            flint_abort(); \
                        } \
                    }

#define TEST_OP1(gr_op, arf_op, opname) \
                    if (gr_op(z, x, ctx) == GR_SUCCESS) \
                    { \
                        nfloat_get_arf(az2, z, ctx); \
                        arf_op(az, ax, prec, rnd); \
                        if (arf_cmp(az2, az) == ((rnd == ARF_RND_FLOOR) ? 1 : -1)) \
                        { \
                            flint_printf("FAIL: %s\n", opname); \
                            flint_printf("prec = %wd, direction = %s\n", prec, (rnd == ARF_RND_FLOOR) ? "floor" : "ceil"); \
                            flint_printf("x = %{gr}\n", x, ctx); \
                            flint_printf("z = %{gr}\n", z, ctx); \
                            flint_printf("ax = "); arf_printd(ax, prec / 3.32 + 1); flint_printf("\n"); \
                            flint_printf("az = "); arf_printd(az, prec / 3.32 + 1); flint_printf("\n"); \
                            flint_printf("az2 = "); arf_printd(az2, prec / 3.32 + 1); flint_printf("\n"); \
                            flint_abort(); \
                        } \
                    }

                    TEST_OP(gr_add, arf_add, "add")
                    TEST_OP(gr_sub, arf_sub, "sub")
                    TEST_OP(gr_mul, arf_mul, "mul")
                    TEST_OP(gr_div, arf_div, "div")
                    TEST_OP(gr_addmul, arf_addmul, "addmul")
                    TEST_OP(gr_submul, arf_submul, "submul")

                    TEST_OP1(gr_sqr, _arf_sqr, "sqr")
                    TEST_OP1(gr_inv, _arf_inv, "inv")
                    TEST_OP1(gr_sqrt, arf_sqrt, "sqrt")
                    TEST_OP1(gr_rsqrt, arf_rsqrt, "rsqrt")
#undef TEST_OP
#undef TEST_OP1

                    gr_ptr X, Y;
                    arf_ptr aX, aY;
                    slong j, N = n_randint(state, 5);
                    int initial = n_randint(state, 2);
                    int subtract = n_randint(state, 2);
                    int reverse = n_randint(state, 2);
                    int status;

                    X = gr_heap_init_vec(N, ctx);
                    Y = gr_heap_init_vec(N, ctx);
                    aX = _arf_vec_init(N);
                    aY = _arf_vec_init(N);

                    for (j = 0; j < N; j++)
                    {
                        GR_MUST_SUCCEED(gr_randtest(GR_ENTRY(X, j, ctx->sizeof_elem), state, ctx));
                        GR_MUST_SUCCEED(gr_randtest(GR_ENTRY(Y, j, ctx->sizeof_elem), state, ctx));
                        nfloat_get_arf(aX + j, GR_ENTRY(X, j, ctx->sizeof_elem), ctx);
                        nfloat_get_arf(aY + j, GR_ENTRY(Y, j, ctx->sizeof_elem), ctx);
                    }

                    nfloat_get_arf(ax, x, ctx);

                    if (reverse)
                        status = _gr_vec_dot_rev(z, initial ? x : NULL, subtract, X, Y, N, ctx);
                    else
                        status = _gr_vec_dot(z, initial ? x : NULL, subtract, X, Y, N, ctx);

                    if (status == GR_SUCCESS)
                    {
                        nfloat_get_arf(az2, z, ctx);

                        if (reverse)
                            arf_dot(az, initial ? ax : NULL, subtract, aX, 1, aY + N - 1, -1, N, prec, rnd);
                        else
                            arf_dot(az, initial ? ax : NULL, subtract, aX, 1, aY, 1, N, prec, rnd);

                        if (arf_cmp(az2, az) == ((rnd == ARF_RND_FLOOR) ? 1 : -1))
                        {
                            flint_printf("FAIL: dot\n");
                            flint_printf("prec = %wd, direction = %s\n", prec, (rnd == ARF_RND_FLOOR) ? "floor" : "ceil");
                            flint_printf("X = %{gr*}\n", X, N, ctx);
                            flint_printf("Y = %{gr*}\n", Y, N, ctx);
                            flint_printf("x = "); arf_printd(x, prec / 3.32 + 1); flint_printf("\n");
                            flint_printf("az = "); arf_printd(az, prec / 3.32 + 1); flint_printf("\n");
                            flint_printf("az2 = "); arf_printd(az2, prec / 3.32 + 1); flint_printf("\n");
                            flint_abort();
                        }
                    }

                    gr_heap_clear_vec(X, N, ctx);
                    gr_heap_clear_vec(Y, N, ctx);
                    _arf_vec_clear(aX, N);
                    _arf_vec_clear(aY, N);
                }

                arf_clear(ax);
                arf_clear(ay);
                arf_clear(az);
                arf_clear(az2);

                gr_heap_clear(x, ctx);
                gr_heap_clear(y, ctx);
                gr_heap_clear(z, ctx);
            }

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

                    gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_re, ctx2, tol, state, 0);
                    gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_im, ctx2, tol, state, 0);
                    gr_test_approx_unary_op(ctx, (gr_method_unary_op) gr_conj, ctx2, tol, state, 0);

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
    }

    TEST_FUNCTION_END(state);
}
