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

int
gr_test_cmp_fun(gr_ctx_t R, gr_method_binary_op_get_int op, gr_ctx_t R_ref, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    int cmp1, cmp2;
    gr_ptr a, b, a_ref, b_ref;

    GR_TMP_INIT2(a, b, R);
    GR_TMP_INIT2(a_ref, b_ref, R_ref);

    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);

    status |= gr_set_other(a_ref, a, R, R_ref);
    status |= gr_set_other(b_ref, b, R, R_ref);

    status |= op(&cmp1, a, b, R);
    status |= op(&cmp2, a_ref, b_ref, R_ref);

    if (status == GR_SUCCESS && cmp1 != cmp2)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        gr_ctx_println(R_ref);
        flint_printf("a = "); gr_println(a, R);
        flint_printf("b = "); gr_println(b, R);
        flint_printf("cmp1 = %d\n", cmp1);
        flint_printf("cmp2 = %d\n", cmp2);
        flint_printf("\n");
    }

    if (status == GR_TEST_FAIL)
        flint_abort();

    GR_TMP_CLEAR2(a, b, R);
    GR_TMP_CLEAR2(a_ref, b_ref, R_ref);

    return status;
}

int
gr_test_approx_unary_op(gr_ctx_t R, gr_method_unary_op op, gr_ctx_t R_ref, gr_srcptr rel_tol, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    int alias;
    int cmp;
    gr_ptr a, b, a_ref, b_ref, rel_err;

    GR_TMP_INIT2(a, b, R);
    GR_TMP_INIT3(a_ref, b_ref, rel_err, R_ref);

    alias = n_randint(state, 2);

    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);

    status |= gr_set_other(a_ref, a, R, R_ref);

    if (status == GR_SUCCESS)
    {
        if (alias == 0)
        {
            status |= op(b, a, R);
            status |= op(b_ref, a_ref, R_ref);
        }
        else
        {
            status |= gr_set(b, a, R);
            status |= op(b, b, R);
            status |= op(b_ref, a_ref, R_ref);
        }

        if (status == GR_SUCCESS)
        {
            status |= gr_set_other(rel_err, b, R, R_ref);

            status |= gr_sub(rel_err, b_ref, rel_err, R_ref);
            status |= gr_div(rel_err, rel_err, b_ref, R_ref);
            status |= gr_abs(rel_err, rel_err, R_ref);
            status |= gr_cmp(&cmp, rel_err, rel_tol, R_ref);

            if (status == GR_SUCCESS && cmp > 0)
                status = GR_TEST_FAIL;

            if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
            {
                flint_printf("\n");
                gr_ctx_println(R);
                gr_ctx_println(R_ref);
                flint_printf("alias: %d\n", alias);
                flint_printf("Computed:\n");
                flint_printf("a = "); gr_println(a, R);
                flint_printf("op(a) = "); gr_println(b, R);
                flint_printf("Reference:\n");
                flint_printf("a = "); gr_println(a_ref, R_ref);
                flint_printf("op(a) = "); gr_println(b_ref, R_ref);
                flint_printf("\nrel_err = "); gr_println(rel_err, R_ref);
                flint_printf("\nrel_tol = "); gr_println(rel_tol, R_ref);
                flint_printf("\n");
            }
        }
    }

    if (status == GR_TEST_FAIL)
        flint_abort();

    GR_TMP_CLEAR2(a, b, R);
    GR_TMP_CLEAR3(a_ref, b_ref, rel_err, R_ref);

    return status;
}

int
gr_test_approx_binary_op(gr_ctx_t R, gr_method_binary_op op, gr_ctx_t R_ref, gr_srcptr rel_tol, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    int alias;
    int cmp;
    gr_ptr a, b, c, a_ref, b_ref, c_ref, rel_err;

    GR_TMP_INIT3(a, b, c, R);
    GR_TMP_INIT4(a_ref, b_ref, c_ref, rel_err, R_ref);

    alias = n_randint(state, 5);

    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_randtest(c, state, R);

    status |= gr_set_other(a_ref, a, R, R_ref);
    status |= gr_set_other(b_ref, b, R, R_ref);

    if (status == GR_SUCCESS)
    {
        if (alias == 0)
        {
            status |= op(c, a, b, R);
            status |= op(c_ref, a_ref, b_ref, R_ref);
        }
        else if (alias == 1)
        {
            status |= gr_set(c, a, R);
            status |= op(c, c, b, R);
            status |= op(c_ref, a_ref, b_ref, R_ref);
        }
        else if (alias == 2)
        {
            status |= gr_set(c, b, R);
            status |= op(c, a, c, R);
            status |= op(c_ref, a_ref, b_ref, R_ref);
        }
        else if (alias == 3)
        {
            status |= gr_set(b, a, R);
            status |= gr_set(b_ref, a_ref, R_ref);
            status |= op(c, a, a, R);
            status |= op(c_ref, a_ref, a_ref, R_ref);
        }
        else
        {
            status |= gr_set(b, a, R);
            status |= gr_set(c, a, R);
            status |= gr_set(b_ref, a_ref, R_ref);
            status |= op(c, c, c, R);
            status |= op(c_ref, a_ref, a_ref, R_ref);
        }

        if (status == GR_SUCCESS)
        {
            status |= gr_set_other(rel_err, c, R, R_ref);

            status |= gr_sub(rel_err, c_ref, rel_err, R_ref);
            status |= gr_div(rel_err, rel_err, c_ref, R_ref);
            status |= gr_abs(rel_err, rel_err, R_ref);
            status |= gr_cmp(&cmp, rel_err, rel_tol, R_ref);

            if (status == GR_SUCCESS && cmp > 0)
                status = GR_TEST_FAIL;

            if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
            {
                flint_printf("\n");
                gr_ctx_println(R);
                gr_ctx_println(R_ref);
                flint_printf("alias: %d\n", alias);
                flint_printf("Computed:\n");
                flint_printf("a = "); gr_println(a, R);
                flint_printf("b = "); gr_println(b, R);
                flint_printf("a (op) b = "); gr_println(c, R);
                flint_printf("Reference:\n");
                flint_printf("a = "); gr_println(a_ref, R_ref);
                flint_printf("b = "); gr_println(b_ref, R_ref);
                flint_printf("a (op) b = "); gr_println(c_ref, R_ref);
                flint_printf("\nrel_err = "); gr_println(rel_err, R_ref);
                flint_printf("\nrel_tol = "); gr_println(rel_tol, R_ref);
                flint_printf("\n");
            }
        }
    }

    if (status == GR_TEST_FAIL)
        flint_abort();

    GR_TMP_CLEAR3(a, b, c, R);
    GR_TMP_CLEAR4(a_ref, b_ref, c_ref, rel_err, R_ref);

    return status;
}

int
gr_test_approx_binary_op_type_variants(gr_ctx_t R,
    const char * opname,
    gr_method_binary_op gr_op,
    gr_method_binary_op_ui gr_op_ui,
    gr_method_binary_op_si gr_op_si,
    gr_method_binary_op_fmpz gr_op_fmpz,
    gr_method_binary_op_fmpq gr_op_fmpq,
    int fused,
    int small_test_values,
    gr_srcptr rel_tol, flint_rand_t state, int test_flags)
{
    int status, alias, which;
    gr_ptr x, y, xy1, xy2, err;
    ulong uy;
    slong sy;
    fmpz_t zy;
    fmpq_t qy;

    GR_TMP_INIT5(x, y, xy1, xy2, err, R);
    fmpz_init(zy);
    fmpq_init(qy);

    which = 0;

    if (small_test_values)
    {
        uy = n_randint(state, 6);
        sy = -5 + (slong) n_randint(state, 11);
        fmpz_randtest(zy, state, 3);
        fmpq_randtest(qy, state, 3);
    }
    else
    {
        uy = n_randtest(state);
        sy = (slong) n_randtest(state);

        if (n_randint(state, 10) == 0)
        {
            fmpz_randtest(zy, state, 10000);
            fmpq_randtest(qy, state, 200);
        }
        else
        {
            fmpz_randtest(zy, state, 100);
            fmpq_randtest(qy, state, 100);
        }
    }

    for (which = 0; which < 4; which++)
    {
        status = GR_SUCCESS;
        alias = n_randint(state, 2);

        GR_MUST_SUCCEED(gr_randtest(x, state, R));
        GR_MUST_SUCCEED(gr_randtest(y, state, R));
        GR_MUST_SUCCEED(gr_randtest(xy1, state, R));

        if (fused && alias)
            GR_MUST_SUCCEED(gr_set(xy2, x, R));
        else if (fused)
            GR_MUST_SUCCEED(gr_set(xy2, xy1, R));
        else
            GR_MUST_SUCCEED(gr_randtest(xy2, state, R));

        if (alias)
            GR_MUST_SUCCEED(gr_set(xy1, x, R));

        if (which == 0)
        {
            if (alias)
                status |= gr_op_ui(xy1, xy1, uy, R);
            else
                status |= gr_op_ui(xy1, x, uy, R);

            status |= gr_set_ui(y, uy, R);

            if (status == GR_SUCCESS)
                status |= gr_op(xy2, x, y, R);
        }
        else if (which == 1)
        {
            if (alias)
                status |= gr_op_si(xy1, xy1, sy, R);
            else
                status |= gr_op_si(xy1, x, sy, R);
            status |= gr_set_si(y, sy, R);

            if (status == GR_SUCCESS)
                status |= gr_op(xy2, x, y, R);
        }
        else if (which == 2)
        {
            if (alias)
                status |= gr_op_fmpz(xy1, xy1, zy, R);
            else
                status |= gr_op_fmpz(xy1, x, zy, R);
            status |= gr_set_fmpz(y, zy, R);

            if (status == GR_SUCCESS)
                status |= gr_op(xy2, x, y, R);
        }
        else
        {
            if (alias)
                status |= gr_op_fmpq(xy1, xy1, qy, R);
            else
                status |= gr_op_fmpq(xy1, x, qy, R);
            status |= gr_set_fmpq(y, qy, R);

            if (status == GR_SUCCESS)
                status |= gr_op(xy2, x, y, R);
        }

        if (status == GR_SUCCESS)
        {
            int cmp;

            status |= gr_sub(err, xy1, xy2, R);
            status |= gr_div(err, err, xy2, R);
            status |= gr_abs(err, err, R);
            status |= gr_cmp(&cmp, err, rel_tol, R);

            if (status == GR_SUCCESS && cmp > 0)
            {
                status = GR_TEST_FAIL;
                break;
            }
        }
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("%s\n", opname);
        gr_ctx_println(R);
        flint_printf("which: %d\n", which);
        flint_printf("alias: %d\n", alias);
        flint_printf("x = "); gr_println(x, R);
        flint_printf("y = "); gr_println(y, R);
        flint_printf("y (op) y (1) = "); gr_println(xy1, R);
        flint_printf("x (op) y (2) = "); gr_println(xy2, R);
        flint_printf("err = "); gr_println(err, R);
        flint_printf("tol = "); gr_println(rel_tol, R);
        flint_printf("\n");
    }

    if (status == GR_TEST_FAIL)
        flint_abort();

    GR_TMP_CLEAR5(x, y, xy1, xy2, err, R);

    fmpz_clear(zy);
    fmpq_clear(qy);

    return status;
}

void
nfloat_test_dot(flint_rand_t state, slong iters, gr_ctx_t ctx)
{
    slong iter, i, len, prec, ebits;
    gr_ptr s0, vec1, vec2, res;
    int subtract, initial, reverse;
    arf_ptr avec1, avec2;
    arf_t as0, ares, ares2, amag, err, t;

    prec = NFLOAT_CTX_PREC(ctx);

    for (iter = 0; iter < iters; iter++)
    {
        len = n_randint(state, 30);

        if (n_randint(state, 2))
            ebits = 1;
        else
            ebits = 10;

        initial = n_randint(state, 2);
        subtract = n_randint(state, 2);
        reverse = n_randint(state, 2);

        vec1 = gr_heap_init_vec(len, ctx);
        vec2 = gr_heap_init_vec(len, ctx);
        s0 = gr_heap_init(ctx);
        res = gr_heap_init(ctx);

        avec1 = _arf_vec_init(len);
        avec2 = _arf_vec_init(len);
        arf_init(as0);
        arf_init(ares);
        arf_init(ares2);
        arf_init(amag);
        arf_init(t);
        arf_init(err);

        if (initial)
        {
            arf_randtest(as0, state, prec, ebits);
            GR_MUST_SUCCEED(nfloat_set_arf(s0, as0, ctx));

            arf_set(ares, as0);
            arf_abs(t, as0);
            arf_add(amag, amag, t, prec, ARF_RND_DOWN);
        }

        for (i = 0; i < len; i++)
        {
            if (i > 0 && n_randint(state, 2))
            {
                arf_set(avec1 + i, avec1 + len - 1 - i);
                arf_neg(avec2 + i, avec2 + len - 1 - i);
            }
            else
            {
                arf_randtest(avec1 + i, state, prec, ebits);
                arf_randtest(avec2 + i, state, prec, ebits);
            }
        }

        for (i = 0; i < len; i++)
        {
            arf_mul(t, avec1 + i, avec2 + (reverse ? len - 1 - i : i), 2 * prec, ARF_RND_DOWN);

            if (subtract)
                arf_sub(ares, ares, t, 2 * prec, ARF_RND_DOWN);
            else
                arf_add(ares, ares, t, 2 * prec, ARF_RND_DOWN);

            arf_abs(t, t);
            arf_add(amag, amag, t, prec, ARF_RND_DOWN);
        }

        /* tolerance */
        arf_mul_2exp_si(t, amag, -prec + 3);

        for (i = 0; i < len; i++)
        {
            GR_MUST_SUCCEED(nfloat_set_arf(GR_ENTRY(vec1, i, ctx->sizeof_elem), avec1 + i, ctx));
            GR_MUST_SUCCEED(nfloat_set_arf(GR_ENTRY(vec2, i, ctx->sizeof_elem), avec2 + i, ctx));
        }

        if (reverse)
            GR_MUST_SUCCEED(_nfloat_vec_dot_rev(res, initial ? s0 : NULL, subtract, vec1, vec2, len, ctx));
        else
            GR_MUST_SUCCEED(_nfloat_vec_dot(res, initial ? s0 : NULL, subtract, vec1, vec2, len, ctx));

        GR_MUST_SUCCEED(nfloat_get_arf(ares2, res, ctx));

        arf_sub(err, ares, ares2, prec, ARF_RND_DOWN);
        arf_abs(err, err);

        if (arf_cmpabs(err, t) > 0)
        {
            flint_printf("FAIL: dot\n");
            gr_ctx_println(ctx);

            flint_printf("reverse = %d, subtract = %d\n", reverse, subtract);

            if (initial)
            {
                flint_printf("\n\ninitial = ");
                arf_printd(as0, 2 + prec / 3.33);
            }

            flint_printf("\n\nvec1 = ");
            _gr_vec_print(vec1, len, ctx);
            flint_printf("\n\nvec2 = ");
            _gr_vec_print(vec2, len, ctx);
            flint_printf("\n\nares = \n");
            arf_printd(ares, 2 + prec / 3.33);
            flint_printf("\n\nares2 = \n");
            arf_printd(ares2, 2 + prec / 3.33);
            flint_printf("\n\ntol = \n");
            arf_printd(t, 10);
            flint_printf("\n\nerr = \n");
            arf_printd(err, 10);
            flint_printf("\n\n");

            flint_abort();
        }

        gr_heap_clear_vec(vec1, len, ctx);
        gr_heap_clear_vec(vec2, len, ctx);
        gr_heap_clear(s0, ctx);
        gr_heap_clear(res, ctx);

        _arf_vec_clear(avec1, len);
        _arf_vec_clear(avec2, len);
        arf_clear(as0);
        arf_clear(ares);
        arf_clear(ares2);
        arf_clear(amag);
        arf_clear(t);
        arf_clear(err);
    }
}

TEST_FUNCTION_START(nfloat, state)
{
    gr_ctx_t ctx;
    gr_ctx_t ctx2;
    slong prec;

    for (prec = NFLOAT_MIN_LIMBS * FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        nfloat_ctx_init(ctx, prec, 0);

        gr_test_floating_point(ctx, 100 * flint_test_multiplier(), 0);

        nfloat_test_dot(state, (prec <= 128 ? 10000 : 100) * flint_test_multiplier(), ctx);

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

            gr_ctx_init_real_arb(ctx2, prec + 32);

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

    TEST_FUNCTION_END(state);
}
