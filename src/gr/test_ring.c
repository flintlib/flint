/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "profiler.h"
#include "long_extras.h"

#include "fexpr.h"
#include "fmpq.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"

typedef int ((*gr_test_function)(gr_ctx_t, flint_rand_t, int));

int
gr_test_binary_op_aliasing(gr_ctx_t R, int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t), flint_rand_t state, int test_flags)
{
    int status, alias;
    gr_ptr x, y, xy1, xy2;

    GR_TMP_INIT4(x, y, xy1, xy2, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status = GR_SUCCESS;
    status |= gr_op(xy1, x, y, R);

    alias = n_randint(state, 4);
    switch (alias)
    {
        case 0:
            status |= gr_set(xy2, x, R);
            status |= gr_op(xy1, x, y, R);
            status |= gr_op(xy2, xy2, y, R);
            break;
        case 1:
            status |= gr_set(xy2, y, R);
            status |= gr_op(xy1, x, y, R);
            status |= gr_op(xy2, x, xy2, R);
            break;
        case 2:
            status |= gr_set(y, x, R);
            status |= gr_op(xy1, x, y, R);
            status |= gr_op(xy2, x, x, R);
            break;
        default:
            status |= gr_set(y, x, R);
            status |= gr_set(xy2, x, R);
            status |= gr_op(xy1, x, y, R);
            status |= gr_op(xy2, xy2, xy2, R);
    }

    if (status == GR_SUCCESS && gr_equal(xy1, xy2, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("alias: %d\n", alias);
        flint_printf("x = "); gr_println(x, R);
        flint_printf("y = "); gr_println(y, R);
        flint_printf("y (op) y (1) = "); gr_println(xy1, R);
        flint_printf("x (op) y (2) = "); gr_println(xy2, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy1, xy2, R);

    return status;
}

int
gr_test_set_ui(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr xa, xb, xc, xa_xb;
    ulong a, b, c;

    do {
        a = n_randtest(state);
        b = n_randtest(state);
        c = a + b;
    } while (c < a);

    GR_TMP_INIT4(xa, xb, xc, xa_xb, R);

    GR_MUST_SUCCEED(gr_randtest(xa, state, R));

    status = GR_SUCCESS;
    status |= gr_set_ui(xa, a, R);
    status |= gr_set_ui(xb, b, R);
    status |= gr_set_ui(xc, c, R);
    status |= gr_add(xa_xb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xc, xa_xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && a == 1 && gr_is_one(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && a == 0 && gr_is_zero(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && b == 1 && gr_is_one(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && b == 0 && gr_is_zero(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("a = %wu\n", a);
        flint_printf("b = %wu\n", b);
        flint_printf("c = %wu\n", c);
        flint_printf("xa = "); gr_println(xa, R);
        flint_printf("xb = "); gr_println(xb, R);
        flint_printf("xc = "); gr_println(xc, R);
        flint_printf("xa + xb = "); gr_println(xa_xb, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(xa, xb, xc, xa_xb, R);

    return status;
}

int
gr_test_set_si(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr xa, xb, xc, xa_xb;
    slong a, b, c;

    do {
        a = z_randtest(state);
        b = z_randtest(state);
    }
    while (z_add_checked(&c, a, b));

    GR_TMP_INIT4(xa, xb, xc, xa_xb, R);

    GR_MUST_SUCCEED(gr_randtest(xa, state, R));

    status = GR_SUCCESS;
    status |= gr_set_si(xa, a, R);
    status |= gr_set_si(xb, b, R);
    status |= gr_set_si(xc, c, R);
    status |= gr_add(xa_xb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xc, xa_xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if (status == GR_SUCCESS && a == 1 && gr_is_one(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && a == 0 && gr_is_zero(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && b == 1 && gr_is_one(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && b == 0 && gr_is_zero(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("a = %wd\n", a);
        flint_printf("b = %wd\n", b);
        flint_printf("c = %wd\n", c);
        flint_printf("xa = "); gr_println(xa, R);
        flint_printf("xb = "); gr_println(xb, R);
        flint_printf("xc = "); gr_println(xc, R);
        flint_printf("xa + xb = "); gr_println(xa_xb, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(xa, xb, xc, xa_xb, R);

    return status;
}

int
gr_test_set_fmpz(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr xa, xb, xc, xa_xb;
    fmpz_t a, b, c;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);

    fmpz_randtest(a, state, 100);
    fmpz_randtest(b, state, 100);
    fmpz_add(c, a, b);

    GR_TMP_INIT4(xa, xb, xc, xa_xb, R);

    GR_MUST_SUCCEED(gr_randtest(xa, state, R));

    status = GR_SUCCESS;
    status |= gr_set_fmpz(xa, a, R);
    status |= gr_set_fmpz(xb, b, R);
    status |= gr_set_fmpz(xc, c, R);
    status |= gr_add(xa_xb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xc, xa_xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if (status == GR_SUCCESS && fmpz_is_one(a) && gr_is_one(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpz_is_zero(a) && gr_is_zero(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpz_is_one(b) && gr_is_one(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpz_is_zero(b) && gr_is_zero(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("a = "); fmpz_print(a); flint_printf("\n");
        flint_printf("b = "); fmpz_print(b); flint_printf("\n");
        flint_printf("c = "); fmpz_print(c); flint_printf("\n");
        flint_printf("xa = "); gr_println(xa, R);
        flint_printf("xb = "); gr_println(xb, R);
        flint_printf("xc = "); gr_println(xc, R);
        flint_printf("xa + xb = "); gr_println(xa_xb, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(xa, xb, xc, xa_xb, R);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);

    return status;
}

int
gr_test_set_fmpq(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr xa, xb, xc, xa_xb;
    fmpq_t a, b, c;

    fmpq_init(a);
    fmpq_init(b);
    fmpq_init(c);

    fmpq_randtest(a, state, 100);
    fmpq_randtest(b, state, 100);
    fmpq_add(c, a, b);

    GR_TMP_INIT4(xa, xb, xc, xa_xb, R);

    GR_MUST_SUCCEED(gr_randtest(xa, state, R));

    status = GR_SUCCESS;
    status |= gr_set_fmpq(xa, a, R);
    status |= gr_set_fmpq(xb, b, R);
    status |= gr_set_fmpq(xc, c, R);
    status |= gr_add(xa_xb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xc, xa_xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if (status == GR_SUCCESS && fmpq_is_one(a) && gr_is_one(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpq_is_zero(a) && gr_is_zero(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpq_is_one(b) && gr_is_one(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpq_is_zero(b) && gr_is_zero(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("a = "); fmpq_print(a); flint_printf("\n");
        flint_printf("b = "); fmpq_print(b); flint_printf("\n");
        flint_printf("c = "); fmpq_print(c); flint_printf("\n");
        flint_printf("xa = "); gr_println(xa, R);
        flint_printf("xb = "); gr_println(xb, R);
        flint_printf("xc = "); gr_println(xc, R);
        flint_printf("xa + xb = "); gr_println(xa_xb, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(xa, xb, xc, xa_xb, R);

    fmpq_clear(a);
    fmpq_clear(b);
    fmpq_clear(c);

    return status;
}

int
gr_test_get_ui(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y;
    ulong a;

    GR_TMP_INIT2(x, y, R);

    status = GR_SUCCESS;

    if (n_randint(state, 2))
    {
        GR_MUST_SUCCEED(gr_randtest(x, state, R));
    }
    else
    {
        a = n_randtest(state);
        status |= gr_set_ui(x, a, R);
        a = n_randtest(state);
    }

    status |= gr_get_ui(&a, x, R);

    if (status == GR_SUCCESS)
    {
        status |= gr_set_ui(y, a, R);

        if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
            status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = "); gr_println(x, R);
        flint_printf("a = %wu\n", a);
        flint_printf("y = "); gr_println(y, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR2(x, y, R);

    return status;
}

int
gr_test_get_si(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y;
    slong a;

    GR_TMP_INIT2(x, y, R);

    status = GR_SUCCESS;

    if (n_randint(state, 2))
    {
        GR_MUST_SUCCEED(gr_randtest(x, state, R));
    }
    else
    {
        a = n_randtest(state);
        status |= gr_set_si(x, a, R);
        a = n_randtest(state);
    }

    status |= gr_get_si(&a, x, R);

    if (status == GR_SUCCESS)
    {
        status |= gr_set_si(y, a, R);

        if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
            status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = "); gr_println(x, R);
        flint_printf("a = %wd\n", a);
        flint_printf("y = "); gr_println(y, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR2(x, y, R);

    return status;
}

int
gr_test_get_fmpz(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y;
    fmpz_t a;

    GR_TMP_INIT2(x, y, R);
    fmpz_init(a);

    status = GR_SUCCESS;

    if (n_randint(state, 2))
    {
        GR_MUST_SUCCEED(gr_randtest(x, state, R));
    }
    else
    {
        fmpz_randtest(a, state, 100);
        status |= gr_set_fmpz(x, a, R);
        fmpz_randtest(a, state, 100);
    }

    status |= gr_get_fmpz(a, x, R);

    if (status == GR_SUCCESS)
    {
        status |= gr_set_fmpz(y, a, R);

        if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
            status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = "); gr_println(x, R);
        flint_printf("a = \n"); fmpz_print(a); flint_printf("\n");
        flint_printf("y = "); gr_println(y, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR2(x, y, R);
    fmpz_clear(a);

    return status;
}

int
gr_test_get_fmpq(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y;
    fmpq_t a;

    GR_TMP_INIT2(x, y, R);
    fmpq_init(a);

    status = GR_SUCCESS;

    if (n_randint(state, 2))
    {
        GR_MUST_SUCCEED(gr_randtest(x, state, R));
    }
    else
    {
        fmpq_randtest(a, state, 100);
        status |= gr_set_fmpq(x, a, R);
        fmpq_randtest(a, state, 100);
    }

    status |= gr_get_fmpq(a, x, R);

    if (status == GR_SUCCESS)
    {
        status |= gr_set_fmpq(y, a, R);

        if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
            status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = "); gr_println(x, R);
        flint_printf("a = \n"); fmpq_print(a); flint_printf("\n");
        flint_printf("y = "); gr_println(y, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR2(x, y, R);
    fmpq_clear(a);

    return status;
}

int
gr_test_get_fmpz_2exp_fmpz(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y;
    fmpz_t a, b;

    GR_TMP_INIT2(x, y, R);
    fmpz_init(a);
    fmpz_init(b);

    status = GR_SUCCESS;

    if (n_randint(state, 2))
    {
        GR_MUST_SUCCEED(gr_randtest(x, state, R));
    }
    else
    {
        fmpz_randtest(a, state, 100);
        fmpz_randtest(b, state, 10);
        status |= gr_set_fmpz_2exp_fmpz(x, a, b, R);
        fmpz_randtest(a, state, 100);
        fmpz_randtest(b, state, 10);
    }

    status |= gr_get_fmpz_2exp_fmpz(a, b, x, R);

    if (status == GR_SUCCESS)
    {
        status |= gr_set_fmpz_2exp_fmpz(y, a, b, R);

        if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
            status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = "); gr_println(x, R);
        flint_printf("a = \n"); fmpz_print(a); flint_printf("\n");
        flint_printf("b = \n"); fmpz_print(b); flint_printf("\n");
        flint_printf("y = "); gr_println(y, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR2(x, y, R);
    fmpz_clear(a);
    fmpz_clear(b);

    return status;
}

int
gr_test_get_set_fexpr(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y;
    fexpr_t expr;

    GR_TMP_INIT2(x, y, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    fexpr_init(expr);

    status |= gr_get_fexpr(expr, x, R);

    if (status == GR_SUCCESS)
    {
        fexpr_vec_t inp;
        gr_vec_t out;

        fexpr_vec_init(inp, 0);
        gr_vec_init(out, 0, R);

        status |= gr_set_fexpr(y, inp, out, expr, R);

        fexpr_vec_clear(inp);
        gr_vec_clear(out, R);

        if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("get_set_fexpr\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        fexpr_print(expr); flint_printf("\n");
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR2(x, y, R);

    fexpr_clear(expr);

    return status;
}

int
gr_test_ctx_get_str(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    char * s;

    status = gr_ctx_get_str(&s, R);

    if (status != GR_SUCCESS)
    {
        status = GR_TEST_FAIL;
        flint_printf("ctx_get_str\n");
    }

    flint_free(s);

    return status;
}

int
gr_test_get_set_str(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y;
    char * s;

    GR_TMP_INIT2(x, y, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status |= gr_get_str(&s, x, R);

    if (status == GR_SUCCESS)
    {
        status |= gr_set_str(y, s, R);

        if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }
    }
    else
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("get_set_str\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        if (s == NULL)
            flint_printf("(NULL)\n");
        else
            flint_printf("%s\n", s);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("\n");
    }

    flint_free(s);

    GR_TMP_CLEAR2(x, y, R);

    return status;
}

int
gr_test_set_other(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, z, xy, x2, y2, z2, t2;
    gr_ctx_t R2;

    gr_ctx_init_random(R2, state);

    GR_TMP_INIT4(x, y, z, xy, R);
    GR_TMP_INIT4(x2, y2, z2, t2, R2);

    status |= gr_randtest(x2, state, R2);
    status |= gr_randtest(y2, state, R2);
    status |= gr_add(z2, x2, y2, R2);

    status |= gr_set_other(x, x2, R2, R);
    status |= gr_set_other(y, y2, R2, R);
    status |= gr_set_other(z, z2, R2, R);

    status |= gr_add(xy, x, y, R);

    /* check that the conversion is homomorphic */
    if (status == GR_SUCCESS && gr_equal(xy, z, R) == T_FALSE)
        status = GR_TEST_FAIL;

    /* test the reverse conversions */
    if (status == GR_SUCCESS && gr_set_other(t2, x, R, R2) == GR_SUCCESS && gr_equal(x2, t2, R2) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && gr_set_other(t2, y, R, R2) == GR_SUCCESS && gr_equal(y2, t2, R2) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && gr_set_other(t2, xy, R, R2) == GR_SUCCESS && gr_equal(z2, t2, R2) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("gr_test_set_other\n");
        gr_ctx_println(R);
        gr_ctx_println(R2);
        flint_printf("x2 = \n"); gr_println(x2, R2);
        flint_printf("y2 = \n"); gr_println(y2, R2);
        flint_printf("z2 = \n"); gr_println(z2, R2);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("z = \n"); gr_println(z, R);
        flint_printf("xy = \n"); gr_println(xy, R);
        flint_printf("t2 = \n"); gr_println(t2, R2);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, z, xy, R);
    GR_TMP_CLEAR4(x2, y2, z2, t2, R2);

    gr_ctx_clear(R2);

    return status;
}

int
gr_test_mul_2exp_si(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, r1, r2;
    slong y;

    GR_TMP_INIT3(x, r1, r2, R);

    status = GR_SUCCESS;

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(r1, state, R));
    y = n_randint(state, 200) - 100;

    if (n_randint(state, 2))
    {
        status |= gr_mul_2exp_si(r1, x, y, R);
    }
    else
    {
        status |= gr_set(r1, x, R);
        status |= gr_mul_2exp_si(r1, r1, y, R);
    }

    if (n_randint(state, 2))
    {
        status |= gr_set_ui(r2, 2, R);
        status |= gr_pow_si(r2, r2, y, R);
        status |= gr_mul(r2, x, r2, R);
    }
    else
    {
        status |= gr_set_ui(r2, 2, R);
        status |= gr_pow_si(r2, r2, -y, R);
        status |= gr_div(r2, x, r2, R);
    }

    if (status == GR_SUCCESS && gr_equal(r1, r2, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = "); gr_println(x, R);
        flint_printf("y = %wd\n", y);
        flint_printf("r1 = "); gr_println(r1, R);
        flint_printf("r2 = "); gr_println(r2, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, r1, r2, R);

    return status;
}

int
gr_test_mul_2exp_fmpz(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, r1, r2;
    fmpz_t y;

    GR_TMP_INIT3(x, r1, r2, R);
    fmpz_init(y);

    status = GR_SUCCESS;

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(r1, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
        fmpz_randbits(y, state, 100);
    else
        fmpz_randbits(y, state, 8);

    if (n_randint(state, 2))
    {
        status |= gr_mul_2exp_fmpz(r1, x, y, R);
    }
    else
    {
        status |= gr_set(r1, x, R);
        status |= gr_mul_2exp_fmpz(r1, r1, y, R);
    }

    if (n_randint(state, 2))
    {
        status |= gr_set_ui(r2, 2, R);
        status |= gr_pow_fmpz(r2, r2, y, R);
        status |= gr_mul(r2, x, r2, R);
    }
    else
    {
        status |= gr_set_ui(r2, 2, R);
        fmpz_neg(y, y);
        status |= gr_pow_fmpz(r2, r2, y, R);
        fmpz_neg(y, y);
        status |= gr_div(r2, x, r2, R);
    }

    if (status == GR_SUCCESS && gr_equal(r1, r2, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = "); gr_println(x, R);
        flint_printf("y = "); fmpz_print(y); flint_printf("\n");
        flint_printf("r1 = "); gr_println(r1, R);
        flint_printf("r2 = "); gr_println(r2, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, r1, r2, R);
    fmpz_clear(y);

    return status;
}

int
gr_test_binary_op_type_variants(gr_ctx_t R,
    const char * opname,
    int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    int (*gr_op_ui)(gr_ptr, gr_srcptr, ulong, gr_ctx_t),
    int (*gr_op_si)(gr_ptr, gr_srcptr, slong, gr_ctx_t),
    int (*gr_op_fmpz)(gr_ptr, gr_srcptr, const fmpz_t, gr_ctx_t),
    int (*gr_op_fmpq)(gr_ptr, gr_srcptr, const fmpq_t, gr_ctx_t),
    int fused,
    int small_test_values,
    flint_rand_t state, int test_flags)
{
    int status, alias, which;
    gr_ptr x, y, xy1, xy2;
    ulong uy;
    slong sy;
    fmpz_t zy;
    fmpq_t qy;

    GR_TMP_INIT4(x, y, xy1, xy2, R);
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
        fmpz_randtest(zy, state, 100);
        fmpq_randtest(qy, state, 100);
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

        if (status == GR_SUCCESS && gr_equal(xy1, xy2, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
            break;
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
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy1, xy2, R);

    fmpz_clear(zy);
    fmpq_clear(qy);

    return status;
}

int
gr_test_binary_op_associative(gr_ctx_t R, int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t), flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, z;
    gr_ptr xy, yz, xy_z, x_yz;

    GR_TMP_INIT3(x, y, z, R);
    GR_TMP_INIT4(xy, yz, xy_z, x_yz, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));
    GR_MUST_SUCCEED(gr_randtest(xy, state, R));
    GR_MUST_SUCCEED(gr_randtest(yz, state, R));

    status = GR_SUCCESS;
    status |= gr_op(xy, x, y, R);
    status |= gr_op(yz, y, z, R);
    status |= gr_op(xy_z, xy, z, R);
    status |= gr_op(x_yz, x, yz, R);

    if (status == GR_SUCCESS && gr_equal(xy_z, x_yz, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("z = \n"); gr_println(z, R);
        flint_printf("x (op) y = \n"); gr_println(xy, R);
        flint_printf("y (op) z = \n"); gr_println(yz, R);
        flint_printf("(x (op) y) (op) z = \n"); gr_println(xy_z, R);
        flint_printf("x (op) (y (op) z) = \n"); gr_println(x_yz, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, y, z, R);
    GR_TMP_CLEAR4(xy, yz, xy_z, x_yz, R);

    return status;
}

int
gr_test_binary_op_commutative(gr_ctx_t R, int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t), flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, xy, yx;

    GR_TMP_INIT4(x, y, xy, yx, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status = GR_SUCCESS;
    status |= gr_op(xy, x, y, R);
    status |= gr_op(yx, y, x, R);

    if (status == GR_SUCCESS && gr_equal(xy, yx, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("y (op) y = \n"); gr_println(xy, R);
        flint_printf("y (op) x = \n"); gr_println(yx, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy, yx, R);

    return status;
}

/*
test x op (y op2 z) = (x op y) op2 (x op z)
*/
int
gr_test_binary_op_left_distributive(gr_ctx_t R,
    int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    int (*gr_op2)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, z, yz, x_yz, xy, xz, xy_xz;

    GR_TMP_INIT4(x, y, z, yz, R);
    GR_TMP_INIT4(x_yz, xy, xz, xy_xz, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));

    status = GR_SUCCESS;
    status |= gr_op2(yz, y, z, R);
    status |= gr_op(x_yz, x, yz, R);
    status |= gr_op(xy, x, y, R);
    status |= gr_op(xz, x, z, R);
    status |= gr_op2(xy_xz, xy, xz, R);

    if (status == GR_SUCCESS && gr_equal(x_yz, xy_xz, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("y (op2) z = \n"); gr_println(yz, R);
        flint_printf("x (op) (y (op2) z) = \n"); gr_println(x_yz, R);
        flint_printf("x (op) y = \n"); gr_println(xy, R);
        flint_printf("x (op) z = \n"); gr_println(xz, R);
        flint_printf("(x op y) (op2) (x op z) = \n"); gr_println(xy_xz, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, z, yz, R);
    GR_TMP_CLEAR4(x_yz, xy, xz, xy_xz, R);

    return status;
}

/*
test (y op2 z) op x = (y op x) op2 (z op x)
*/
int
gr_test_binary_op_right_distributive(gr_ctx_t R,
    int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    int (*gr_op2)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, z, yz, yz_x, yx, zx, yx_zx;

    GR_TMP_INIT4(x, y, z, yz, R);
    GR_TMP_INIT4(yz_x, yx, zx, yx_zx, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));

    status = GR_SUCCESS;
    status |= gr_op2(yz, y, z, R);
    status |= gr_op(yz_x, yz, x, R);
    status |= gr_op(yx, y, x, R);
    status |= gr_op(zx, z, x, R);
    status |= gr_op2(yx_zx, yx, zx, R);

    if (status == GR_SUCCESS && gr_equal(yz_x, yx_zx, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("y (op2) z = \n"); gr_println(yz, R);
        flint_printf("(y (op2) z) op x = \n"); gr_println(yz_x, R);
        flint_printf("y (op) x = \n"); gr_println(yz, R);
        flint_printf("z (op) x = \n"); gr_println(zx, R);
        flint_printf("(y op x) (op2) (z op x) = \n"); gr_println(yx_zx, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, z, yz, R);
    GR_TMP_CLEAR4(yz_x, yx, zx, yx_zx, R);

    return status;
}

int
gr_test_init_clear(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, b, c, d, e;

    status = GR_SUCCESS;

    GR_TMP_INIT(a, R);
    status |= gr_randtest(a, state, R);
    GR_TMP_CLEAR(a, R);

    GR_TMP_INIT2(a, b, R);
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    GR_TMP_CLEAR2(a, b, R);

    GR_TMP_INIT3(a, b, c, R);
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_randtest(c, state, R);
    GR_TMP_CLEAR3(a, b, c, R);

    GR_TMP_INIT4(a, b, c, d, R);
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_randtest(c, state, R);
    status |= gr_randtest(d, state, R);
    GR_TMP_CLEAR4(a, b, c, d, R);

    GR_TMP_INIT5(a, b, c, d, e, R);
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_randtest(c, state, R);
    status |= gr_randtest(d, state, R);
    status |= gr_randtest(e, state, R);
    GR_TMP_CLEAR5(a, b, c, d, e, R);

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    return status;
}

int
gr_test_equal(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, b;
    truth_t equal0, equal1;

    status = GR_SUCCESS;

    GR_TMP_INIT2(a, b, R);

    status |= gr_randtest(a, state, R);
    status |= gr_set(b, a, R);

    equal0 = gr_equal(a, a, R);
    equal1 = gr_equal(a, b, R);

    if (status == GR_SUCCESS && (equal0 == T_FALSE || equal1 == T_FALSE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if (status == GR_TEST_FAIL)
    {
        flint_printf("FAIL: equal\n");
        gr_ctx_println(R);
        flint_printf("a = "); gr_println(a, R);
        flint_printf("(a == a) = "); truth_println(equal0);
        flint_printf("b = "); gr_println(b, R);
        flint_printf("(a == b) = "); truth_println(equal1);
    }

    GR_TMP_CLEAR2(a, b, R);

    return status;
}

int
gr_test_swap(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, b, c, d;
    truth_t equal0, equal1, equal2, equal3, equal4;

    status = GR_SUCCESS;

    GR_TMP_INIT4(a, b, c, d, R);

    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_set(c, a, R);
    status |= gr_set(d, b, R);
    gr_swap(a, a, R);

    equal0 = gr_equal(a, c, R);

    gr_swap(a, b, R);
    equal1 = gr_equal(b, c, R);
    equal2 = gr_equal(a, d, R);

    gr_swap(a, b, R);
    equal3 = gr_equal(a, c, R);
    equal4 = gr_equal(b, d, R);

    if (status == GR_SUCCESS &&
        (equal0 == T_FALSE || equal1 == T_FALSE ||
         equal2 == T_FALSE || equal3 == T_FALSE ||
         equal4 == T_FALSE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    GR_TMP_CLEAR4(a, b, c, d, R);

    return status;
}

int
gr_test_zero_one(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a;
    truth_t equal;

    status = GR_SUCCESS;

    GR_TMP_INIT(a, R);

    status |= gr_randtest(a, state, R);
    status |= gr_zero(a, R);
    equal = gr_is_zero(a, R);
    if (status == GR_SUCCESS && equal == T_FALSE)
        status = GR_TEST_FAIL;

    status |= gr_randtest(a, state, R);
    status |= gr_one(a, R);
    equal = gr_is_one(a, R);
    if (status == GR_SUCCESS && equal == T_FALSE)
    {
        flint_printf("FAILL is_one\n");
        gr_ctx_println(R);
        gr_println(a, R);
        status = GR_TEST_FAIL;
    }

    status |= gr_randtest(a, state, R);
    status |= gr_one(a, R);
    status |= gr_neg(a, a, R);
    equal = gr_is_neg_one(a, R);
    if (status == GR_SUCCESS && equal == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    GR_TMP_CLEAR(a, R);

    return status;
}

int
gr_test_randtest_not_zero(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a;
    truth_t is_zero;

    status = GR_SUCCESS;

    GR_TMP_INIT(a, R);

    status |= gr_randtest_not_zero(a, state, R);

    if (status == GR_SUCCESS)
    {
        is_zero = gr_is_zero(a, R);

        if (status == GR_SUCCESS && is_zero == T_TRUE)
            status = GR_TEST_FAIL;
    }

    GR_TMP_CLEAR(a, R);

    return status;
}


int
gr_test_one(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a;
    truth_t equal;

    status = GR_SUCCESS;

    GR_TMP_INIT(a, R);

    status |= gr_randtest(a, state, R);
    status |= gr_one(a, R);
    equal = gr_is_one(a, R);
    if (status == GR_SUCCESS && equal == T_FALSE)
        status = GR_TEST_FAIL;

    status |= gr_randtest(a, state, R);
    status |= gr_one(a, R);
    status |= gr_inv(a, a, R);
    if (status == GR_SUCCESS && equal == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    GR_TMP_CLEAR(a, R);

    return status;
}

int
gr_test_add_associative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_associative(R, gr_add, state, test_flags);
}

int
gr_test_neg(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, xy;

    GR_TMP_INIT3(x, y, xy, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(xy, state, R));

    status = GR_SUCCESS;

    /* check x + (-x) = 0 */
    status |= gr_neg(y, x, R);
    status |= gr_add(xy, x, y, R);

    if (status == GR_SUCCESS && gr_is_zero(xy, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("x + y = \n"); gr_println(xy, R);
        flint_printf("\n");
    }

    /* check -(-x) = x, plus aliasing */
    status |= gr_neg(y, y, R);

    if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, y, xy, R);

    return status;
}

int
gr_test_add_commutative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_commutative(R, gr_add, state, test_flags);
}

int
gr_test_add_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_aliasing(R, gr_add, state, test_flags);
}

int
gr_test_add_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R, "add",
        gr_add, gr_add_ui, gr_add_si, gr_add_fmpz, gr_add_fmpq,
            0, 0, state, test_flags);
}

int
gr_test_sub_equal_neg_add(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, neg_y, x_sub_y, x_neg_y;

    GR_TMP_INIT5(x, y, neg_y, x_sub_y, x_neg_y, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(neg_y, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_sub_y, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_neg_y, state, R));

    status = GR_SUCCESS;
    status |= gr_sub(x_sub_y, x, y, R);
    status |= gr_neg(neg_y, y, R);
    status |= gr_add(x_neg_y, x, neg_y, R);

    if (status == GR_SUCCESS && gr_equal(x_sub_y, x_neg_y, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("-y = \n"); gr_println(neg_y, R);
        flint_printf("x - y = \n"); gr_println(x_sub_y, R);
        flint_printf("x + (-y) = \n"); gr_println(x_neg_y, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR5(x, y, neg_y, x_sub_y, x_neg_y, R);

    return status;
}

int
gr_test_sub_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_aliasing(R, gr_sub, state, test_flags);
}

int
gr_test_sub_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R, "sub",
        gr_sub, gr_sub_ui, gr_sub_si, gr_sub_fmpz, gr_sub_fmpq,
            0, 0, state, test_flags);
}

int
gr_test_mul_associative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_associative(R, gr_mul, state, test_flags);
}

int
gr_test_mul_commutative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_commutative(R, gr_mul, state, test_flags);
}

int
gr_test_mul_left_distributive(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_left_distributive(R, gr_mul, gr_add, state, test_flags);
}

int
gr_test_mul_right_distributive(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_right_distributive(R, gr_mul, gr_add, state, test_flags);
}


int
gr_test_mul_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_aliasing(R, gr_mul, state, test_flags);
}

int
gr_test_mul_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R, "mul",
        gr_mul, gr_mul_ui, gr_mul_si, gr_mul_fmpz, gr_mul_fmpq,
            0, 0, state, test_flags);
}

int
gr_test_addmul_submul(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    int which;
    gr_ptr x, y, z, t;

    GR_TMP_INIT4(x, y, z, t, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));

    status = GR_SUCCESS;

    which = n_randint(state, 10);

    switch (which)
    {
        case 0:
            status |= gr_mul(t, y, z, R);
            status |= gr_add(t, x, t, R);
            status |= gr_addmul(x, y, z, R);
            break;
        case 1:
            status |= gr_mul(t, y, y, R);
            status |= gr_add(t, x, t, R);
            status |= gr_addmul(x, y, y, R);
            break;
        case 2:
            status |= gr_mul(t, x, z, R);
            status |= gr_add(t, x, t, R);
            status |= gr_addmul(x, x, z, R);
            break;
        case 3:
            status |= gr_mul(t, y, x, R);
            status |= gr_add(t, x, t, R);
            status |= gr_addmul(x, y, x, R);
            break;
        case 4:
            status |= gr_mul(t, x, x, R);
            status |= gr_add(t, x, t, R);
            status |= gr_addmul(x, x, x, R);
            break;
        case 5:
            status |= gr_mul(t, y, z, R);
            status |= gr_sub(t, x, t, R);
            status |= gr_submul(x, y, z, R);
            break;
        case 6:
            status |= gr_mul(t, y, y, R);
            status |= gr_sub(t, x, t, R);
            status |= gr_submul(x, y, y, R);
            break;
        case 7:
            status |= gr_mul(t, x, z, R);
            status |= gr_sub(t, x, t, R);
            status |= gr_submul(x, x, z, R);
            break;
        case 8:
            status |= gr_mul(t, y, x, R);
            status |= gr_sub(t, x, t, R);
            status |= gr_submul(x, y, x, R);
            break;
        case 9:
            status |= gr_mul(t, x, x, R);
            status |= gr_sub(t, x, t, R);
            status |= gr_submul(x, x, x, R);
            break;
        default:
            break;
    }

    if (status == GR_SUCCESS && gr_equal(x, t, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("which = %d\n", which);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("z = \n"); gr_println(z, R);
        flint_printf("t = \n"); gr_println(t, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, z, t, R);

    return status;
}

int
gr_test_addmul_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R, "addmul",
        gr_addmul, gr_addmul_ui, gr_addmul_si, gr_addmul_fmpz, gr_addmul_fmpq,
            1, 0, state, test_flags);
}

int
gr_test_submul_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R, "submul",
        gr_submul, gr_submul_ui, gr_submul_si, gr_submul_fmpz, gr_submul_fmpq,
            1, 0, state, test_flags);
}

int
gr_test_div_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_aliasing(R, gr_div, state, test_flags);
}

int
gr_test_div_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R, "div",
        gr_div, gr_div_ui, gr_div_si, gr_div_fmpz, gr_div_fmpq,
            0, 0, state, test_flags);
}

int
gr_test_is_invertible(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    truth_t invertible = T_UNKNOWN;
    gr_ptr x, x_inv;

    GR_TMP_INIT2(x, x_inv, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));

    status = GR_SUCCESS;
    status |= gr_inv(x_inv, x, R);

    if (status != GR_UNABLE)
    {
        invertible = gr_is_invertible(x, R);

        if ((status == GR_SUCCESS && invertible == T_FALSE) ||
            (status == GR_DOMAIN && invertible == T_TRUE))
        {
            status = GR_TEST_FAIL;
        }
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("x ^ -1 = \n"); gr_println(x_inv, R);
        flint_printf("status = %d, invertible = %d\n", status, invertible);
        flint_printf("\n");
    }

    GR_TMP_CLEAR2(x, x_inv, R);

    return status;
}

int
gr_test_inv_involution(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, x_inv, x_inv_inv;

    GR_TMP_INIT3(x, x_inv, x_inv_inv, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_inv, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_inv_inv, state, R));

    status = GR_SUCCESS;
    status |= gr_inv(x_inv, x, R);
    status |= gr_inv(x_inv_inv, x_inv, R);

    if (status == GR_SUCCESS && gr_equal(x, x_inv_inv, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("x ^ -1 = \n"); gr_println(x_inv, R);
        flint_printf("(x ^ -1) ^ -1 = \n"); gr_println(x_inv_inv, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, x_inv, x_inv_inv, R);

    return status;
}

int
gr_test_inv_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    truth_t equal1, equal2;
    gr_ptr x, x_inv, x_inv_x, x_x_inv;

    GR_TMP_INIT4(x, x_inv, x_inv_x, x_x_inv, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_inv, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_inv_x, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_x_inv, state, R));

    /* todo: split status */
    status = GR_SUCCESS;
    status |= gr_inv(x_inv, x, R);
    status |= gr_mul(x_inv_x, x_inv, x, R);
    status |= gr_mul(x_x_inv, x, x_inv, R);
    equal1 = gr_is_one(x_inv_x, R);
    equal2 = gr_is_one(x_x_inv, R);

    if (status == GR_SUCCESS && (equal1 == T_FALSE || equal2 == T_FALSE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("x ^ -1 = \n"); gr_println(x_inv, R);
        flint_printf("(x ^ -1) * x = \n"); gr_println(x_inv_x, R);
        flint_printf("x * (x ^ -1) = \n"); gr_println(x_x_inv, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, x_inv, x_inv_x, x_x_inv, R);

    return status;
}

int
gr_test_div_right_distributive(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_right_distributive(R, gr_div, gr_add, state, test_flags);
}

int
gr_test_div_then_mul(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, xy, xyy;

    GR_TMP_INIT4(x, y, xy, xyy, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(xy, state, R));
    GR_MUST_SUCCEED(gr_randtest(xyy, state, R));

    status = GR_SUCCESS;
    status |= gr_div(xy, x, y, R);
    status |= gr_mul(xyy, xy, y, R);

    if (status == GR_SUCCESS && gr_equal(x, xyy, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("x / y = \n"); gr_println(xy, R);
        flint_printf("(x / y) * y = \n"); gr_println(xyy, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy, xyy, R);

    return status;
}

int
gr_test_mul_then_div(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, xy, xyy;

    GR_TMP_INIT4(x, y, xy, xyy, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(xy, state, R));
    GR_MUST_SUCCEED(gr_randtest(xyy, state, R));

    status = GR_SUCCESS;
    status |= gr_mul(xy, x, y, R);
    status |= gr_div(xyy, xy, y, R);

    if (status == GR_SUCCESS && gr_equal(x, xyy, R) == T_FALSE && gr_ctx_is_integral_domain(R) == T_TRUE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("x * y = \n"); gr_println(xy, R);
        flint_printf("(x * y) / y = \n"); gr_println(xyy, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy, xyy, R);

    return status;
}

int
gr_test_divexact(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    int aliasing;
    gr_ptr x, y, xy, q;

    GR_TMP_INIT4(x, y, xy, q, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    aliasing = n_randint(state, 3);

    status = GR_SUCCESS;
    status |= gr_mul(xy, x, y, R);

    if (aliasing == 0)
    {
        status |= gr_divexact(q, xy, y, R);
    }
    else if (aliasing == 1)
    {
        status |= gr_set(q, xy, R);
        status |= gr_divexact(q, q, y, R);
    }
    else
    {
        status |= gr_set(q, y, R);
        status |= gr_divexact(q, xy, q, R);
    }

    if (status == GR_SUCCESS && gr_equal(q, x, R) == T_FALSE && gr_ctx_is_integral_domain(R) == T_TRUE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("aliasing = %d\n", aliasing);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("xy = \n"); gr_println(xy, R);
        flint_printf("q = \n"); gr_println(q, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy, q, R);

    return status;
}

/* todo: fmpq */
int
gr_test_divexact_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status, alias, which;
    gr_ptr x, xy, q;
    ulong uy;
    slong sy;
    fmpz_t zy;

    GR_TMP_INIT3(x, xy, q, R);
    fmpz_init(zy);

    which = 0;

    uy = n_randtest(state);
    sy = (slong) n_randtest(state);
    fmpz_randtest(zy, state, 100);

    for (which = 0; which < 4; which++)
    {
        GR_MUST_SUCCEED(gr_randtest(x, state, R));
        GR_MUST_SUCCEED(gr_randtest(q, state, R));

        status = GR_SUCCESS;
        alias = n_randint(state, 2);

        if (which == 0)
        {
            status |= gr_mul_ui(xy, x, uy, R);

            if (alias)
            {
                status |= gr_set(q, xy, R);
                status |= gr_divexact_ui(q, q, uy, R);
            }
            else
            {
                status |= gr_divexact_ui(q, xy, uy, R);
            }
        }
        else if (which == 1)
        {
            status |= gr_mul_si(xy, x, sy, R);

            if (alias)
            {
                status |= gr_set(q, xy, R);
                status |= gr_divexact_si(q, q, sy, R);
            }
            else
            {
                status |= gr_divexact_si(q, xy, sy, R);
            }
        }
        else
        {
            status |= gr_mul_fmpz(xy, x, zy, R);

            if (alias)
            {
                status |= gr_set(q, xy, R);
                status |= gr_divexact_fmpz(q, q, zy, R);
            }
            else
            {
                status |= gr_divexact_fmpz(q, xy, zy, R);
            }
        }

        if (status == GR_SUCCESS && gr_equal(q, x, R) == T_FALSE && gr_ctx_is_integral_domain(R) == T_TRUE)
        {
            status = GR_TEST_FAIL;
            break;
        }
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("which: %d\n", which);
        flint_printf("alias: %d\n", alias);
        flint_printf("x = "); gr_println(x, R);
        flint_printf("xy = "); gr_println(xy, R);
        flint_printf("q = "); gr_println(q, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, xy, q, R);

    fmpz_clear(zy);

    return status;
}

int
gr_test_div_nonunique(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    int status2 = GR_SUCCESS;
    int status3 = GR_SUCCESS;
    int status4 = GR_SUCCESS;
    int status5 = GR_SUCCESS;
    gr_ptr x, y, xy, z, q;

    GR_TMP_INIT5(x, y, xy, z, q, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status |= gr_mul(xy, x, y, R);

    if (status == GR_SUCCESS)
    {
        status2 = gr_div_nonunique(q, xy, x, R);

        if (status2 == GR_DOMAIN)
        {
            status = GR_TEST_FAIL;
        }
        else if (status2 == GR_SUCCESS)
        {
            status2 = gr_mul(z, q, x, R);
            if (status2 == GR_SUCCESS && gr_equal(z, xy, R) == T_FALSE)
                status = GR_TEST_FAIL;
        }

        status3 = gr_div_nonunique(q, xy, y, R);

        if (status3 == GR_DOMAIN)
        {
            status = GR_TEST_FAIL;
        }
        else if (status3 == GR_SUCCESS)
        {
            status3 = gr_mul(z, q, y, R);
            if (status3 == GR_SUCCESS && gr_equal(z, xy, R) == T_FALSE)
                status = GR_TEST_FAIL;
        }

        status4 = gr_div_nonunique(z, x, y, R);

        if (status4 == GR_DOMAIN)
        {
            /* if claimed non-divisible, verify that div claims the same */
            status5 = gr_div(z, x, y, R);
            if (status5 == GR_SUCCESS)
                status = GR_TEST_FAIL;
        }
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("div_nonunique\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("xy = \n"); gr_println(xy, R);
        flint_printf("z = \n"); gr_println(z, R);
        flint_printf("status = %d, %d, %d, %d, %d\n", status,
            status2, status3, status4, status5);
        flint_printf("\n");
    }

    GR_TMP_CLEAR5(x, y, xy, z, q, R);

    return status;
}

int
gr_test_div_nonunique_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_aliasing(R, gr_div_nonunique, state, test_flags);
}

int
gr_test_divides(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    int status2 = GR_SUCCESS;
    int status3 = GR_SUCCESS;
    int status4 = GR_SUCCESS;
    truth_t x_divides = T_UNKNOWN, y_divides = T_UNKNOWN;
    gr_ptr x, y, xy, z;

    GR_TMP_INIT4(x, y, xy, z, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status |= gr_mul(xy, x, y, R);

    if (status == GR_SUCCESS)
    {
        x_divides = gr_divides(x, xy, R);
        y_divides = gr_divides(y, xy, R);

        if (x_divides == T_FALSE || y_divides == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }

        if (gr_ctx_is_integral_domain(R) == T_TRUE)
        {
            if (gr_is_zero(x, R) == T_FALSE)
            {
                status2 = gr_divexact(z, xy, x, R);
                if (status2 == GR_DOMAIN)
                    status = GR_TEST_FAIL;
            }

            if (gr_is_zero(y, R) == T_FALSE)
            {
                status3 = gr_divexact(z, xy, y, R);
                if (status3 == GR_DOMAIN)
                    status = GR_TEST_FAIL;
            }
        }
    }

    if (status == GR_SUCCESS)
    {
        truth_t d = gr_divides(x, y, R);

        /* if claimed non-divisible, verify that div claims the same */
        if (d == T_FALSE)
        {
            status4 = gr_div(z, y, x, R);

            if (status4 == GR_SUCCESS)
                status = GR_TEST_FAIL;
        }
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("divides\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("xy = \n"); gr_println(xy, R);
        flint_printf("x divides = "); truth_println(x_divides);
        flint_printf("y divides = "); truth_println(y_divides);
        flint_printf("status = %d, %d, %d, %d\n", status,
            status2, status3, status4);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy, z, R);

    return status;
}

int
gr_test_pow_ui_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    ulong a, b;
    gr_ptr x, xa, xb, xab, xaxb;

    GR_TMP_INIT5(x, xa, xb, xab, xaxb, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(xb, state, R));
    GR_MUST_SUCCEED(gr_randtest(xab, state, R));
    GR_MUST_SUCCEED(gr_randtest(xaxb, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
    {
        do {
            a = n_randtest(state);
            b = n_randtest(state);
        } while (a + b < a);
    }
    else
    {
        a = n_randtest(state) % 20;
        b = n_randtest(state) % 20;
    }

    status = GR_SUCCESS;

    status |= gr_pow_ui(xa, x, a, R);
    status |= gr_pow_ui(xb, x, b, R);
    status |= gr_pow_ui(xab, x, a + b, R);
    status |= gr_mul(xaxb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xab, xaxb, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("a = %wu\n", a);
        flint_printf("b = %wu\n", b);
        flint_printf("x ^ a = \n"); gr_println(xa, R);
        flint_printf("x ^ b = \n"); gr_println(xb, R);
        flint_printf("x ^ (a + b) = \n"); gr_println(xab, R);
        flint_printf("(x ^ a) * (x ^ b) = \n"); gr_println(xaxb, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR5(x, xa, xb, xab, xaxb, R);

    return status;
}

int
gr_test_pow_ui_base_scalar_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    ulong a;
    slong y;
    gr_ptr x, xa, ya, xya, xaya;

    GR_TMP_INIT3(x, xa, ya, R);
    GR_TMP_INIT2(xya, xaya, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(ya, state, R));

    y = n_randtest(state);

    if (gr_ctx_is_finite(R) == T_TRUE)
        a = n_randtest(state);
    else
        a = n_randtest(state) % 20;

    status = GR_SUCCESS;
    status |= gr_pow_ui(xa, x, a, R);
    status |= gr_set_si(ya, y, R);
    status |= gr_pow_ui(ya, ya, a, R);
    status |= gr_set_si(xya, y, R);
    status |= gr_mul(xya, x, xya, R);   /* todo mul_si */
    status |= gr_pow_ui(xya, xya, a, R);
    status |= gr_mul(xaya, xa, ya, R);

    if (status == GR_SUCCESS && gr_equal(xya, xaya, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = %wd\n", y);
        flint_printf("a = %wu\n", a);
        flint_printf("x ^ a = \n"); gr_println(xa, R);
        flint_printf("y ^ a = \n"); gr_println(ya, R);
        flint_printf("(x * y) ^ a = \n"); gr_println(xya, R);
        flint_printf("(x ^ a) * (y ^ a) = \n"); gr_println(xaya, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, xa, ya, R);
    GR_TMP_CLEAR2(xya, xaya, R);

    return status;
}

int
gr_test_pow_ui_base_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    ulong a;
    gr_ptr x, y, xa, ya, xya, xaya;

    GR_TMP_INIT4(x, y, xa, ya, R);
    GR_TMP_INIT2(xya, xaya, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(ya, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
        a = n_randtest(state);
    else
        a = n_randtest(state) % 20;

    status = GR_SUCCESS;
    status |= gr_pow_ui(xa, x, a, R);
    status |= gr_pow_ui(ya, y, a, R);
    status |= gr_mul(xya, x, y, R);
    status |= gr_pow_ui(xya, xya, a, R);
    status |= gr_mul(xaya, xa, ya, R);

    if (status == GR_SUCCESS && gr_equal(xya, xaya, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("a = %wu\n", a);
        flint_printf("x ^ a = \n"); gr_println(xa, R);
        flint_printf("y ^ a = \n"); gr_println(ya, R);
        flint_printf("(x * y) ^ a = \n"); gr_println(xya, R);
        flint_printf("(x ^ a) * (y ^ a) = \n"); gr_println(xaya, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xa, ya, R);
    GR_TMP_CLEAR2(xya, xaya, R);

    return status;
}

int
gr_test_pow_ui_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    ulong a;
    gr_ptr x, xa1, xa2;

    GR_TMP_INIT3(x, xa1, xa2, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa1, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
        a = n_randtest(state);
    else
        a = n_randtest(state) % 20;

    status = GR_SUCCESS;
    status |= gr_pow_ui(xa1, x, a, R);
    status |= gr_set(xa2, x, R);
    status |= gr_pow_ui(xa2, xa2, a, R);

    if (status == GR_SUCCESS && gr_equal(xa1, xa2, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("a = %wu\n", a);
        flint_printf("x ^ a (1) = \n"); gr_println(xa1, R);
        flint_printf("x ^ a (2) = \n"); gr_println(xa2, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, xa1, xa2, R);

    return status;
}

int
gr_test_pow_fmpz_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    fmpz_t a, b, ab;
    gr_ptr x, xa, xb, xab, xaxb;

    GR_TMP_INIT5(x, xa, xb, xab, xaxb, R);

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(ab);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(xb, state, R));
    GR_MUST_SUCCEED(gr_randtest(xab, state, R));
    GR_MUST_SUCCEED(gr_randtest(xaxb, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
    {
        fmpz_randtest(a, state, 100);
        fmpz_randtest(b, state, 100);
    }
    else if (n_randint(state, 20) == 0)
    {
        if (gr_set_si(x, -1 + (slong) n_randint(state, 3), R) != GR_SUCCESS)
            /* allow using for groups */
            GR_MUST_SUCCEED(gr_one(x, R));
        fmpz_randtest(a, state, 100);
        fmpz_randtest(b, state, 100);
    }
    else
    {
        fmpz_randtest(a, state, 4);
        fmpz_randtest(b, state, 4);
    }

    fmpz_add(ab, a, b);

    status = GR_SUCCESS;

    status |= gr_pow_fmpz(xa, x, a, R);
    status |= gr_pow_fmpz(xb, x, b, R);
    status |= gr_pow_fmpz(xab, x, ab, R);
    status |= gr_mul(xaxb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xab, xaxb, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("a = "); fmpz_print(a); flint_printf("\n");
        flint_printf("b = "); fmpz_print(b); flint_printf("\n");
        flint_printf("x ^ a = \n"); gr_println(xa, R);
        flint_printf("x ^ b = \n"); gr_println(xb, R);
        flint_printf("x ^ (a + b) = \n"); gr_println(xab, R);
        flint_printf("(x ^ a) * (x ^ b) = \n"); gr_println(xaxb, R);
        flint_printf("\n");
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(ab);

    GR_TMP_CLEAR5(x, xa, xb, xab, xaxb, R);

    return status;
}

int
gr_test_pow_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, a, xa1, xa2;

    GR_TMP_INIT4(x, a, xa1, xa2, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa1, state, R));
    GR_MUST_SUCCEED(gr_randtest_small(a, state, R));

    status = GR_SUCCESS;

    switch (n_randint(state, 3))
    {
        case 0:
            status |= gr_set(xa2, x, R);
            status |= gr_pow(xa2, xa2, a, R);
            break;
        case 1:
            status |= gr_set(xa2, a, R);
            status |= gr_pow(xa2, x, xa2, R);
            break;
        default:
            status |= gr_set(xa2, a, R);
            status |= gr_set(x, a, R);
            status |= gr_pow(xa2, xa2, xa2, R);
    }

    status |= gr_pow(xa1, x, a, R);

    if (status == GR_SUCCESS && gr_equal(xa1, xa2, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("a = \n"); gr_println(a, R);
        flint_printf("x ^ a (1) = \n"); gr_println(xa1, R);
        flint_printf("x ^ a (2) = \n"); gr_println(xa2, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, a, xa1, xa2, R);

    return status;
}

int
gr_test_pow_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, b, ab;
    gr_ptr x, xa, xb, xab, xaxb;

    GR_TMP_INIT5(x, xa, xb, xab, xaxb, R);
    GR_TMP_INIT3(a, b, ab, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(xb, state, R));
    GR_MUST_SUCCEED(gr_randtest(xab, state, R));
    GR_MUST_SUCCEED(gr_randtest(xaxb, state, R));

    GR_MUST_SUCCEED(gr_randtest_small(a, state, R));
    GR_MUST_SUCCEED(gr_randtest_small(b, state, R));

    status = GR_SUCCESS;

    status |= gr_add(ab, a, b, R);
    status |= gr_pow(xa, x, a, R);
    status |= gr_pow(xb, x, b, R);
    status |= gr_pow(xab, x, ab, R);
    status |= gr_mul(xaxb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xab, xaxb, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("a = \n"); gr_println(a, R);
        flint_printf("b = \n"); gr_println(b, R);
        flint_printf("a + b = \n"); gr_println(ab, R);
        flint_printf("x ^ a = \n"); gr_println(xa, R);
        flint_printf("x ^ b = \n"); gr_println(xb, R);
        flint_printf("x ^ (a + b) = \n"); gr_println(xab, R);
        flint_printf("(x ^ a) * (x ^ b) = \n"); gr_println(xaxb, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR5(x, xa, xb, xab, xaxb, R);
    GR_TMP_CLEAR3(a, b, ab, R);

    return status;
}

int
gr_test_pow_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R, "pow",
        gr_pow, gr_pow_ui, gr_pow_si, gr_pow_fmpz, gr_pow_fmpq,
            0, 1, state, test_flags);
}

int
gr_test_sqrt(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, y2;
    int perfect;

    GR_TMP_INIT3(x, y, y2, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    perfect = n_randint(state, 2);

    if (perfect)
        status |= gr_sqr(x, x, R);

    if (n_randint(state, 2))
    {
        status |= gr_set(y, x, R);
        status |= gr_sqrt(y, y, R);
    }
    else
    {
        status |= gr_sqrt(y, x, R);
    }

    status |= gr_sqr(y2, y, R);

    if (status == GR_SUCCESS && gr_equal(y2, x, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if (status == GR_DOMAIN && perfect)
    {
        status = GR_TEST_FAIL;
    }

    if (status == GR_SUCCESS && perfect && gr_is_square(x, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("FAIL: sqrt\n");
        flint_printf("R = "); gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("y ^ 2 = \n"); gr_println(y2, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, y, y2, R);

    return status;
}

int
gr_test_rsqrt(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, z;
    int perfect;

    GR_TMP_INIT3(x, y, z, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    perfect = n_randint(state, 2);

    if (perfect)
        status |= gr_sqr(x, x, R);

    if (n_randint(state, 2))
    {
        status |= gr_set(y, x, R);
        status |= gr_rsqrt(y, y, R);
    }
    else
    {
        status |= gr_rsqrt(y, x, R);
    }

    status |= gr_inv(z, y, R);
    status |= gr_sqr(z, z, R);

    if (status == GR_SUCCESS && gr_equal(z, x, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("FAIL: rsqrt\n");
        flint_printf("R = "); gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("1 / y ^ 2 = \n"); gr_println(z, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(x, y, z, R);

    return status;
}

int
gr_test_ordered_ring_cmp(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, z, xz, yz, zero, xy;
    int cmp1, cmp2, cmp3;

    GR_TMP_INIT5(x, y, z, xz, yz, R);
    GR_TMP_INIT2(zero, xy, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));

    /* cmp(x, y) = -cmp(y, x) */
    status |= gr_cmp(&cmp1, x, y, R);
    status |= gr_cmp(&cmp2, y, x, R);

    if (status == GR_SUCCESS && cmp1 != -cmp2)
    {
        status = GR_TEST_FAIL;
    }

    /* x <= y  -->  x + z <= y + z */
    status |= gr_add(xz, x, z, R);
    status |= gr_add(yz, y, z, R);
    status |= gr_cmp(&cmp1, x, y, R);
    status |= gr_cmp(&cmp2, xz, yz, R);

    if (status == GR_SUCCESS && cmp1 != cmp2)
    {
        status = GR_TEST_FAIL;
    }

    /* 0 <= x and 0 <= y --> 0 <= xy */
    status |= gr_cmp(&cmp1, zero, x, R);
    status |= gr_cmp(&cmp2, zero, y, R);
    status |= gr_mul(xy, x, y, R);
    status |= gr_cmp(&cmp3, zero, xy, R);

    if (status == GR_SUCCESS && cmp1 <= 0 && cmp2 <= 0 && cmp3 > 0)
    {
        status = GR_TEST_FAIL;
    }

    if (status & GR_DOMAIN && !(status & GR_UNABLE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("FAIL: ordered_ring_cmp\n");
        flint_printf("R = "); gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("z = \n"); gr_println(z, R);
        flint_printf("x + z = \n"); gr_println(xz, R);
        flint_printf("y + z = \n"); gr_println(yz, R);
        flint_printf("xy = \n"); gr_println(xy, R);
        flint_printf("cmp = %d, %d, %d\n", cmp1, cmp2, cmp3);
        flint_printf("\n");
    }

    GR_TMP_CLEAR5(x, y, z, xz, yz, R);
    GR_TMP_CLEAR2(zero, xy, R);

    return status;
}

int
gr_test_ordered_ring_cmpabs(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, ax, ay;
    int cmp1, cmp2;

    GR_TMP_INIT4(x, y, ax, ay, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status |= gr_abs(ax, x, R);
    status |= gr_abs(ay, y, R);
    status |= gr_cmpabs(&cmp1, x, y, R);
    status |= gr_cmp(&cmp2, ax, ay, R);

    if (status == GR_SUCCESS && cmp1 != cmp2)
    {
        status = GR_TEST_FAIL;
    }

    if (status & GR_DOMAIN && !(status & GR_UNABLE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("R = "); gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("ax = \n"); gr_println(ax, R);
        flint_printf("ay = \n"); gr_println(ay, R);
        flint_printf("cmp = %d, %d\n", cmp1, cmp2);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(x, y, ax, ay, R);

    return status;
}

int
gr_test_complex_parts(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, a, b, ab, i;

    GR_TMP_INIT5(x, a, b, ab, i, R);

    status = gr_i(i, R);

    if (status == GR_SUCCESS)
    {
        int which_test = n_randint(state, 3);

        GR_MUST_SUCCEED(gr_randtest(x, state, R));

        if (which_test == 0)
        {
            /* check x == re(x) + im(x)*i */
            status |= gr_re(a, x, R);
            status |= gr_im(b, x, R);
            status |= gr_mul(ab, b, i, R);
            status |= gr_add(ab, a, ab, R);
        }
        else if (which_test == 1)
        {
            /* check x == abs(x) * sgn(x) */
            status |= gr_abs(a, x, R);
            status |= gr_sgn(b, x, R);
            status |= gr_mul(ab, a, b, R);
        }
        else
        {
            /* check x == re(conj(x)) - im(conj(x))*i */
            status |= gr_conj(a, x, R);
            status |= gr_re(a, a, R);
            status |= gr_conj(b, x, R);
            status |= gr_im(b, b, R);
            status |= gr_mul(ab, b, i, R);
            status |= gr_sub(ab, a, ab, R);
        }

        if (status == GR_SUCCESS && gr_equal(x, ab, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        gr_ctx_println(R);
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("a = \n"); gr_println(a, R);
        flint_printf("b = \n"); gr_println(b, R);
        flint_printf("ab = \n"); gr_println(ab, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR5(x, a, b, ab, i, R);

    return status;
}


int
gr_test_gcd(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, f, g, d, t;
    int aliasing;

    GR_TMP_INIT5(a, f, g, d, t, R);

    status = GR_SUCCESS;
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(f, state, R);
    status |= gr_randtest(g, state, R);
    status |= gr_mul(f, a, f, R);
    status |= gr_mul(g, g, a, R);

    aliasing = n_randint(state, 3);

    if (status == GR_SUCCESS)
    {
        if (aliasing == 0)
        {
            status |= gr_gcd(d, f, g, R);
        }
        else if (aliasing == 1)
        {
            status |= gr_set(d, f, R);
            status |= gr_gcd(d, d, g, R);
        }
        else if (aliasing == 2)
        {
            status |= gr_set(d, g, R);
            status |= gr_gcd(d, f, d, R);
        }

        if (status == GR_SUCCESS && gr_divides(a, d, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("gcd\n");
        gr_ctx_println(R);
        flint_printf("aliasing = %d\n", aliasing);
        flint_printf("a = "); gr_println(a, R);
        flint_printf("f = "); gr_println(f, R);
        flint_printf("g = "); gr_println(g, R);
        flint_printf("gcd = "); gr_println(d, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR5(a, f, g, d, t, R);

    return status;
}

/* verify that LCM(a, b) GCD(a, b) ~ a b */
int
gr_test_lcm(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, b, x, y, ab, xy;
    int aliasing;

    GR_TMP_INIT3(a, b, ab, R);
    GR_TMP_INIT3(x, y, xy, R);

    status = GR_SUCCESS;
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_randtest(x, state, R);
    status |= gr_randtest(y, state, R);

    aliasing = n_randint(state, 3);

    if (aliasing == 0)
    {
        status |= gr_lcm(x, a, b, R);
    }
    else if (aliasing == 1)
    {
        status |= gr_set(x, a, R);
        status |= gr_lcm(x, x, b, R);
    }
    else if (aliasing == 2)
    {
        status |= gr_set(x, b, R);
        status |= gr_lcm(x, a, x, R);
    }

    status |= gr_gcd(y, a, b, R);

    status |= gr_mul(ab, a, b, R);
    status |= gr_mul(xy, x, y, R);

    if (status == GR_SUCCESS)
    {
        if (gr_divides(xy, ab, R) == T_FALSE || gr_divides(ab, xy, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("lcm\n");
        gr_ctx_println(R);
        flint_printf("aliasing = %d\n", aliasing);
        flint_printf("a = "); gr_println(a, R);
        flint_printf("b = "); gr_println(b, R);
        flint_printf("x = "); gr_println(x, R);
        flint_printf("y = "); gr_println(y, R);
        flint_printf("ab = "); gr_println(ab, R);
        flint_printf("xy = "); gr_println(xy, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR3(a, b, ab, R);
    GR_TMP_CLEAR3(x, y, xy, R);

    return status;
}

int
gr_test_numerator_denominator(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, p, q, aq;

    GR_TMP_INIT4(a, p, q, aq, R);

    status = GR_SUCCESS;
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(p, state, R);
    status |= gr_randtest(q, state, R);
    status |= gr_numerator(p, a, R);
    status |= gr_denominator(q, a, R);

    if (status == GR_SUCCESS)
    {
        status |= gr_mul(aq, a, q, R);

        if (status == GR_SUCCESS && gr_equal(aq, p, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("numerator_denominator\n");
        gr_ctx_println(R);
        flint_printf("a = "); gr_println(a, R);
        flint_printf("p = "); gr_println(p, R);
        flint_printf("q = "); gr_println(q, R);
        flint_printf("aq = "); gr_println(aq, R);
        flint_printf("\n");
    }

    GR_TMP_CLEAR4(a, p, q, aq, R);

    return status;
}


int
gr_test_factor(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, c, t, u;
    gr_ctx_t ZZ;
    gr_vec_t fac, exp;
    slong i;

    GR_TMP_INIT4(x, c, t, u, R);
    gr_ctx_init_fmpz(ZZ);

    gr_vec_init(fac, n_randint(state, 3), R);
    gr_vec_init(exp, n_randint(state, 3), ZZ);

    status = GR_SUCCESS;
    status |= gr_randtest_small(x, state, R);

    if (n_randint(state, 2))
    {
        status |= gr_randtest_small(t, state, R);
        status |= gr_mul(x, x, t, R);
    }

    status |= gr_factor(c, fac, exp, x, 0, R);

    if (status == GR_SUCCESS)
    {
        if (fac->length != exp->length)
        {
            status = GR_TEST_FAIL;
        }
        else
        {
            status |= gr_set(u, c, R);

            for (i = 0; i < fac->length; i++)
            {
                status |= gr_pow_fmpz(t, gr_vec_entry_srcptr(fac, i, R), gr_vec_entry_srcptr(exp, i, ZZ), R);
                status |= gr_mul(u, u, t, R);
            }

            if (status == GR_SUCCESS && gr_equal(x, u, R) == T_FALSE)
            {
                status = GR_TEST_FAIL;
            }
        }

        if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
        {
            flint_printf("factor\n");
            flint_printf("x = "); gr_println(x, R);
            flint_printf("c = "); gr_println(c, R);
            flint_printf("fac = "); gr_vec_print(fac, R); flint_printf("\n");
            flint_printf("exp = "); gr_vec_print(exp, ZZ); flint_printf("\n");
            flint_printf("\n");
        }
    }

    GR_TMP_CLEAR4(x, c, t, u, R);
    gr_ctx_clear(ZZ);

    gr_vec_clear(fac, R);
    gr_vec_clear(exp, ZZ);

    return status;
}

int
gr_test_vec_binary_op(gr_ctx_t R, const char * opname, int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    int (*_gr_vec_op)(gr_ptr, gr_srcptr, gr_srcptr, slong, gr_ctx_t), flint_rand_t state, int test_flags)
{
    int status, aliasing;
    slong i, len;
    gr_ptr x, y, xy1, xy2;

    len = n_randint(state, 5);

    GR_TMP_INIT_VEC(x, len, R);
    GR_TMP_INIT_VEC(y, len, R);
    GR_TMP_INIT_VEC(xy1, len, R);
    GR_TMP_INIT_VEC(xy2, len, R);

    GR_MUST_SUCCEED(_gr_vec_randtest(x, state, len, R));

    GR_MUST_SUCCEED(_gr_vec_randtest(y, state, len, R));
    GR_MUST_SUCCEED(_gr_vec_randtest(xy1, state, len, R));
    GR_MUST_SUCCEED(_gr_vec_randtest(xy2, state, len, R));

    status = GR_SUCCESS;

    aliasing = n_randint(state, 4);

    switch (aliasing)
    {
        case 0:
            status |= _gr_vec_set(xy1, x, len, R);
            status |= _gr_vec_op(xy1, xy1, y, len, R);
            break;
        case 1:
            status |= _gr_vec_set(xy1, y, len, R);
            status |= _gr_vec_op(xy1, x, xy1, len, R);
            break;
        case 2:
            status |= _gr_vec_set(y, x, len, R);
            status |= _gr_vec_set(xy1, x, len, R);
            status |= _gr_vec_op(xy1, xy1, xy1, len, R);
            break;
        default:
            status |= _gr_vec_op(xy1, x, y, len, R);
    }

    for (i = 0; i < len; i++)
        status |= gr_op(GR_ENTRY(xy2, i, R->sizeof_elem),
                         GR_ENTRY(x, i, R->sizeof_elem),
                         GR_ENTRY(y, i, R->sizeof_elem), R);

    if (status == GR_SUCCESS && _gr_vec_equal(xy1, xy2, len, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("%s\n", opname);
        gr_ctx_println(R);
        flint_printf("aliasing: %d\n", aliasing);
        _gr_vec_print(xy1, len, R); flint_printf("\n");
        _gr_vec_print(xy2, len, R); flint_printf("\n");
    }

    GR_TMP_CLEAR_VEC(x, len, R);
    GR_TMP_CLEAR_VEC(y, len, R);
    GR_TMP_CLEAR_VEC(xy1, len, R);
    GR_TMP_CLEAR_VEC(xy2, len, R);

    return status;
}

int gr_test_vec_add(gr_ctx_t R, flint_rand_t state, int test_flags) { return gr_test_vec_binary_op(R, "vec_add", gr_add, _gr_vec_add, state, test_flags); }
int gr_test_vec_sub(gr_ctx_t R, flint_rand_t state, int test_flags) { return gr_test_vec_binary_op(R, "vec_sub", gr_sub, _gr_vec_sub, state, test_flags); }
int gr_test_vec_mul(gr_ctx_t R, flint_rand_t state, int test_flags) { return gr_test_vec_binary_op(R, "vec_mul", gr_mul, _gr_vec_mul, state, test_flags); }
int gr_test_vec_div(gr_ctx_t R, flint_rand_t state, int test_flags) { return gr_test_vec_binary_op(R, "vec_div", gr_div, _gr_vec_div, state, test_flags); }
int gr_test_vec_divexact(gr_ctx_t R, flint_rand_t state, int test_flags) { return gr_test_vec_binary_op(R, "vec_divexact", gr_divexact, _gr_vec_divexact, state, test_flags); }
int gr_test_vec_pow(gr_ctx_t R, flint_rand_t state, int test_flags) { return gr_test_vec_binary_op(R, "vec_pow", gr_pow, _gr_vec_pow, state, test_flags); }

int gr_generic_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx);

int
gr_test_vec_dot(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    slong len;
    gr_ptr x, y, a, s, t;
    int initial, alias, subtract, reverse;

    len = n_randint(state, 5);

    initial = n_randint(state, 2);
    alias = n_randint(state, 2);
    subtract = n_randint(state, 2);
    reverse = n_randint(state, 2);

    GR_TMP_INIT_VEC(x, len, R);
    GR_TMP_INIT_VEC(y, len, R);
    GR_TMP_INIT3(a, s, t, R);

    GR_MUST_SUCCEED(_gr_vec_randtest(x, state, len, R));
    GR_MUST_SUCCEED(_gr_vec_randtest(y, state, len, R));
    GR_MUST_SUCCEED(gr_randtest(a, state, R));
    GR_MUST_SUCCEED(gr_randtest(s, state, R));
    GR_MUST_SUCCEED(gr_randtest(t, state, R));

    status = GR_SUCCESS;

    if (initial && alias)
    {
        GR_MUST_SUCCEED(gr_set(s, a, R));
        GR_MUST_SUCCEED(gr_set(t, a, R));
    }

    if (reverse)
    {
        status |= _gr_vec_dot_rev(s, initial ? (alias ? s : a) : NULL, subtract, x, y, len, R);
        status |= _gr_poly_reverse(y, y, len, len, R);
        status |= gr_generic_vec_dot(t, initial ? (alias ? t : a) : NULL, subtract, x, y, len, R);
        status |= _gr_poly_reverse(y, y, len, len, R);
    }
    else
    {
        status |= _gr_vec_dot(s, initial ? (alias ? s : a) : NULL, subtract, x, y, len, R);
        status |= gr_generic_vec_dot(t, initial ? (alias ? t : a) : NULL, subtract, x, y, len, R);
    }

    if (status == GR_SUCCESS && gr_equal(s, t, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("dot\n");
        gr_ctx_println(R);
        flint_printf("alias: %d\n", alias);
        flint_printf("initial: %d\n", initial);
        flint_printf("subtract: %d\n", subtract);
        flint_printf("reverse: %d\n", reverse);
        _gr_vec_print(x, len, R); flint_printf("\n");
        _gr_vec_print(y, len, R); flint_printf("\n");
        gr_println(a, R);
        gr_println(s, R);
        gr_println(t, R);
    }

    GR_TMP_CLEAR_VEC(x, len, R);
    GR_TMP_CLEAR_VEC(y, len, R);
    GR_TMP_CLEAR3(a, s, t, R);

    return status;
}

/* (AB)C = A(BC) */
int
gr_test_mat_mul_classical_associative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_mat_t A, B, C, AB, BC, AB_C, A_BC;
    slong m, n, p, q;

    if (gr_ctx_is_finite(R) == T_TRUE)
    {
        m = n_randint(state, 5);
        n = n_randint(state, 5);
        p = n_randint(state, 5);
        q = n_randint(state, 5);
    }
    else
    {
        m = n_randint(state, 3);
        n = n_randint(state, 3);
        p = n_randint(state, 3);
        q = n_randint(state, 3);
    }

    gr_mat_init(A, m, n, R);
    gr_mat_init(B, n, p, R);
    gr_mat_init(C, p, q, R);
    gr_mat_init(AB, m, p, R);
    gr_mat_init(BC, n, q, R);
    gr_mat_init(AB_C, m, q, R);
    gr_mat_init(A_BC, m, q, R);

    GR_MUST_SUCCEED(gr_mat_randtest(A, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(B, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(C, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(AB, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(BC, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(AB_C, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(A_BC, state, R));

    status = GR_SUCCESS;
    status |= gr_mat_mul_classical(AB, A, B, R);
    status |= gr_mat_mul_classical(BC, B, C, R);
    status |= gr_mat_mul_classical(AB_C, AB, C, R);
    status |= gr_mat_mul_classical(A_BC, A, BC, R);

    if (status == GR_SUCCESS && gr_mat_equal(AB_C, A_BC, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        /* todo: vec print */
        flint_printf("\n");
        flint_printf("A = \n"); gr_mat_print(A, R); flint_printf("\n");
        flint_printf("B = \n"); gr_mat_print(B, R); flint_printf("\n");
        flint_printf("C = \n"); gr_mat_print(C, R); flint_printf("\n");
        flint_printf("AB = \n"); gr_mat_print(AB, R); flint_printf("\n");
        flint_printf("BC = \n"); gr_mat_print(BC, R); flint_printf("\n");
        flint_printf("AB * C = \n"); gr_mat_print(AB_C, R); flint_printf("\n");
        flint_printf("A * BC = \n"); gr_mat_print(A_BC, R); flint_printf("\n");
        flint_printf("\n");
    }

    gr_mat_clear(A, R);
    gr_mat_clear(B, R);
    gr_mat_clear(C, R);
    gr_mat_clear(AB, R);
    gr_mat_clear(BC, R);
    gr_mat_clear(A_BC, R);
    gr_mat_clear(AB_C, R);

    return status;
}

int
gr_test_integral_domain(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, z;

    GR_TMP_INIT3(x, y, z, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    status |= gr_mul(z, x, y, R);

    if (status == GR_SUCCESS && gr_is_zero(x, R) == T_FALSE && gr_is_zero(y, R) == T_FALSE && gr_is_zero(z, R) == T_TRUE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        flint_printf("\n");
        flint_printf("x = \n"); gr_println(x, R);
        flint_printf("y = \n"); gr_println(y, R);
        flint_printf("z = \n"); gr_println(z, R);
        flint_printf("\n");
    }

    if (gr_ctx_is_commutative_ring(R) == T_FALSE)
    {
        flint_printf("integral domain is not a commutative ring\n");
        flint_printf("\n");
        status = GR_TEST_FAIL;
    }

    GR_TMP_CLEAR3(x, y, z, R);

    return status;
}

int
gr_test_field(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, z;

    GR_TMP_INIT3(x, y, z, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));

    if (gr_is_zero(x, R) == T_FALSE)
    {
        if (gr_is_invertible(x, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }

        if (gr_inv(y, x, R) == GR_DOMAIN)
        {
            status = GR_TEST_FAIL;
        }

        if (gr_div(z, y, x, R) == GR_DOMAIN)
        {
            status = GR_TEST_FAIL;
        }

        if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
        {
            flint_printf("\n");
            flint_printf("x = \n"); gr_println(x, R);
            flint_printf("y = \n"); gr_println(y, R);
            flint_printf("z = \n"); gr_println(z, R);
            flint_printf("\n");
        }
    }

    if (gr_ctx_is_commutative_ring(R) == T_FALSE)
    {
        flint_printf("field is not a commutative ring\n");
        flint_printf("\n");
        status = GR_TEST_FAIL;
    }

    if (gr_ctx_is_integral_domain(R) == T_FALSE)
    {
        flint_printf("field is not an integral domain\n");
        flint_printf("\n");
        status = GR_TEST_FAIL;
    }

    GR_TMP_CLEAR3(x, y, z, R);

    return status;
}

void
gr_test_iter(gr_ctx_t R, flint_rand_t state, const char * descr, gr_test_function func, slong iters, int test_flags)
{
    slong iter, count_success, count_unable, count_domain;
    int status;
    timeit_t timer;

    count_success = 0;
    count_unable = 0;
    count_domain = 0;

    if (test_flags & GR_TEST_VERBOSE)
    {
        flint_printf("%s ... ", descr);
        fflush(stdout);
    }

    timeit_start(timer);

    for (iter = 0; iter < iters; iter++)
    {
        /* flint_printf("iter %ld\n", iter); */
        status = func(R, state, test_flags & ~GR_TEST_VERBOSE);

        if (status == GR_SUCCESS)
            count_success++;

        if (status & GR_UNABLE)
            count_unable++;

        if (status & GR_DOMAIN)
            count_domain++;

        if (status & GR_TEST_FAIL)
        {
            flint_throw(FLINT_ERROR, "\nFAIL\n");
        }
    }

    timeit_stop(timer);

    if (test_flags & GR_TEST_VERBOSE)
    {
        flint_printf("PASS   (%wd successful, %wd domain, %wd unable, 0 wrong, %.3g cpu, %.3g wall)\n",
            count_success, count_domain, count_unable, timer->cpu*0.001, timer->wall*0.001);
    }
}

void
gr_test_ring(gr_ctx_t R, slong iters, int test_flags)
{
    timeit_t timer;
    flint_rand_t state;
    slong vec_iters = iters / 10 + 1;

    /* test_flags |= GR_TEST_VERBOSE; */

    if (test_flags & GR_TEST_VERBOSE)
    {
        timeit_start(timer);

        flint_printf("===============================================================================\n");
        flint_printf("Testing "); gr_ctx_println(R);
        flint_printf("-------------------------------------------------------------------------------\n");
    }

    flint_randinit(state);

    /* if (gr_ctx_is_ring(R) != T_TRUE)
        flint_abort(); */

    gr_test_iter(R, state, "ctx_get_str", gr_test_ctx_get_str, 1, test_flags);

    gr_test_iter(R, state, "init/clear", gr_test_init_clear, iters, test_flags);
    gr_test_iter(R, state, "equal", gr_test_equal, iters, test_flags);
    gr_test_iter(R, state, "swap", gr_test_swap, iters, test_flags);
    gr_test_iter(R, state, "zero_one", gr_test_zero_one, iters, test_flags);
    gr_test_iter(R, state, "randtest_not_zero", gr_test_randtest_not_zero, iters, test_flags);
    gr_test_iter(R, state, "neg", gr_test_neg, iters, test_flags);

    gr_test_iter(R, state, "set_ui", gr_test_set_ui, iters, test_flags);
    gr_test_iter(R, state, "set_si", gr_test_set_si, iters, test_flags);
    gr_test_iter(R, state, "set_fmpz", gr_test_set_fmpz, iters, test_flags);
    gr_test_iter(R, state, "set_fmpq", gr_test_set_fmpq, iters, test_flags);
    gr_test_iter(R, state, "set_other", gr_test_set_other, iters, test_flags);

    gr_test_iter(R, state, "get_ui", gr_test_get_ui, iters, test_flags);
    gr_test_iter(R, state, "get_si", gr_test_get_si, iters, test_flags);
    gr_test_iter(R, state, "get_fmpz", gr_test_get_fmpz, iters, test_flags);
    gr_test_iter(R, state, "get_fmpq", gr_test_get_fmpq, iters, test_flags);
    gr_test_iter(R, state, "get_fmpz_2exp_fmpz", gr_test_get_fmpz_2exp_fmpz, iters, test_flags);

    gr_test_iter(R, state, "get_set_fexpr", gr_test_get_set_fexpr, iters, test_flags);
    gr_test_iter(R, state, "get_set_str", gr_test_get_set_str, iters, test_flags);

    gr_test_iter(R, state, "add: associative", gr_test_add_associative, iters, test_flags);
    gr_test_iter(R, state, "add: commutative", gr_test_add_commutative, iters, test_flags);
    gr_test_iter(R, state, "add: aliasing", gr_test_add_aliasing, iters, test_flags);
    gr_test_iter(R, state, "sub: equal neg add", gr_test_sub_equal_neg_add, iters, test_flags);
    gr_test_iter(R, state, "sub: aliasing", gr_test_sub_aliasing, iters, test_flags);

    gr_test_iter(R, state, "add: ui/si/fmpz/fmpq", gr_test_add_type_variants, iters, test_flags);
    gr_test_iter(R, state, "sub: ui/si/fmpz/fmpq", gr_test_sub_type_variants, iters, test_flags);
    gr_test_iter(R, state, "mul: ui/si/fmpz/fmpq", gr_test_mul_type_variants, iters, test_flags);
    gr_test_iter(R, state, "div: ui/si/fmpz/fmpq", gr_test_div_type_variants, iters, test_flags);

    gr_test_iter(R, state, "mul: associative", gr_test_mul_associative, iters, test_flags);
    if (gr_ctx_is_commutative_ring(R) == T_TRUE)
        gr_test_iter(R, state, "mul: commutative", gr_test_mul_commutative, iters, test_flags);
    gr_test_iter(R, state, "mul: aliasing", gr_test_mul_aliasing, iters, test_flags);
    gr_test_iter(R, state, "mul: left distributive", gr_test_mul_left_distributive, iters, test_flags);
    gr_test_iter(R, state, "mul: right distributive", gr_test_mul_right_distributive, iters, test_flags);

    gr_test_iter(R, state, "mul_2exp_si", gr_test_mul_2exp_si, iters, test_flags);
    gr_test_iter(R, state, "mul_2exp_fmpz", gr_test_mul_2exp_fmpz, iters, test_flags);

    gr_test_iter(R, state, "addmul/submul", gr_test_addmul_submul, iters, test_flags);
    gr_test_iter(R, state, "addmul: ui/si/fmpz/fmpq", gr_test_addmul_type_variants, iters, test_flags);
    gr_test_iter(R, state, "submul: ui/si/fmpz/fmpq", gr_test_submul_type_variants, iters, test_flags);

    if (gr_ctx_is_integral_domain(R) == T_TRUE)
        gr_test_iter(R, state, "integral_domain", gr_test_integral_domain, iters, test_flags);

    if (gr_ctx_is_field(R) == T_TRUE)
        gr_test_iter(R, state, "field", gr_test_integral_domain, iters, test_flags);

    if (gr_ctx_is_integral_domain(R) == T_TRUE)
        gr_test_iter(R, state, "div: distributive", gr_test_div_right_distributive, iters, test_flags);

    gr_test_iter(R, state, "div: aliasing", gr_test_div_aliasing, iters, test_flags);

    gr_test_iter(R, state, "div: div then mul", gr_test_div_then_mul, iters, test_flags);
    gr_test_iter(R, state, "div: mul then div", gr_test_mul_then_div, iters, test_flags);

    gr_test_iter(R, state, "div_nonunique", gr_test_div_nonunique, iters, test_flags);
    gr_test_iter(R, state, "div_nonunique: aliasing", gr_test_div_nonunique_aliasing, iters, test_flags);
    gr_test_iter(R, state, "divides", gr_test_divides, iters, test_flags);

    gr_test_iter(R, state, "inv: multiplication", gr_test_inv_multiplication, iters, test_flags);
    gr_test_iter(R, state, "inv: involution", gr_test_inv_involution, iters, test_flags);
    gr_test_iter(R, state, "is_invertible", gr_test_is_invertible, iters, test_flags);

    gr_test_iter(R, state, "divexact", gr_test_divexact, iters, test_flags);
    gr_test_iter(R, state, "divexact: ui/si/fmpz", gr_test_divexact_type_variants, iters, test_flags);

    gr_test_iter(R, state, "pow_ui: exponent addition", gr_test_pow_ui_exponent_addition, iters, test_flags);
    gr_test_iter(R, state, "pow_ui: base scalar multiplication", gr_test_pow_ui_base_scalar_multiplication, iters, test_flags);

    if (gr_ctx_is_commutative_ring(R) == T_TRUE)
        gr_test_iter(R, state, "pow_ui: base multiplication", gr_test_pow_ui_base_multiplication, iters, test_flags);

    gr_test_iter(R, state, "pow_ui: aliasing", gr_test_pow_ui_exponent_addition, iters, test_flags);
    gr_test_iter(R, state, "pow_fmpz: exponent addition", gr_test_pow_fmpz_exponent_addition, iters, test_flags);

    gr_test_iter(R, state, "sqrt", gr_test_sqrt, iters, test_flags & (~GR_TEST_ALWAYS_ABLE));
    gr_test_iter(R, state, "rsqrt", gr_test_rsqrt, iters, test_flags & (~GR_TEST_ALWAYS_ABLE));

    gr_test_iter(R, state, "pow: aliasing", gr_test_pow_aliasing, iters, test_flags & (~GR_TEST_ALWAYS_ABLE));
    gr_test_iter(R, state, "pow: exponent addition", gr_test_pow_exponent_addition, iters, test_flags & (~GR_TEST_ALWAYS_ABLE));
    gr_test_iter(R, state, "pow: ui/si/fmpz/fmpq", gr_test_pow_type_variants, iters, test_flags & (~GR_TEST_ALWAYS_ABLE));

    if (gr_ctx_is_ordered_ring(R) == T_TRUE)
    {
        gr_test_iter(R, state, "ordered_ring_cmp", gr_test_ordered_ring_cmp, iters, test_flags);
        gr_test_iter(R, state, "ordered_ring_cmpabs", gr_test_ordered_ring_cmpabs, iters, test_flags);
    }

    gr_test_iter(R, state, "numerator_denominator", gr_test_numerator_denominator, iters, test_flags);
    gr_test_iter(R, state, "complex_parts", gr_test_complex_parts, iters, test_flags);

    if (gr_ctx_is_unique_factorization_domain(R) == T_TRUE)
    {
        gr_test_iter(R, state, "gcd", gr_test_gcd, iters, test_flags);
        gr_test_iter(R, state, "lcm", gr_test_lcm, iters, test_flags);
        gr_test_iter(R, state, "factor", gr_test_factor, iters, test_flags);
    }

    gr_test_iter(R, state, "vec_add", gr_test_vec_add, vec_iters, test_flags);
    gr_test_iter(R, state, "vec_sub", gr_test_vec_sub, vec_iters, test_flags);
    gr_test_iter(R, state, "vec_mul", gr_test_vec_mul, vec_iters, test_flags);
    gr_test_iter(R, state, "vec_div", gr_test_vec_div, vec_iters, test_flags);
    gr_test_iter(R, state, "vec_divexact", gr_test_vec_divexact, vec_iters, test_flags);
    /* gr_test_iter(R, state, "vec_pow", gr_test_vec_pow, vec_iters, test_flags & (~GR_TEST_ALWAYS_ABLE)); large elements */

    gr_test_iter(R, state, "vec_dot", gr_test_vec_dot, iters, test_flags);

    gr_test_iter(R, state, "mat_mul_classical: associative", gr_test_mat_mul_classical_associative, iters, test_flags);

    flint_randclear(state);

    if (test_flags & GR_TEST_VERBOSE)
    {
        timeit_stop(timer);

        flint_printf("-------------------------------------------------------------------------------\n");
        flint_printf("Tests finished in %.3g cpu, %.3g wall\n", timer->cpu*0.001, timer->wall*0.001);
        flint_printf("===============================================================================\n\n");
    }
}

void
gr_test_multiplicative_group(gr_ctx_t R, slong iters, int test_flags)
{
    timeit_t timer;
    flint_rand_t state;

    /* test_flags |= GR_TEST_VERBOSE; */

    if (test_flags & GR_TEST_VERBOSE)
    {
        timeit_start(timer);

        flint_printf("===============================================================================\n");
        flint_printf("Testing "); gr_ctx_println(R);
        flint_printf("-------------------------------------------------------------------------------\n");
    }

    flint_randinit(state);

    /* if (gr_ctx_is_multiplicative_group(R) != T_TRUE)
        flint_abort(); */

    gr_test_iter(R, state, "ctx_get_str", gr_test_ctx_get_str, 1, test_flags);

    gr_test_iter(R, state, "init/clear", gr_test_init_clear, iters, test_flags);
    gr_test_iter(R, state, "swap", gr_test_swap, iters, test_flags);

    gr_test_iter(R, state, "get_set_fexpr", gr_test_get_set_fexpr, iters, test_flags);
    gr_test_iter(R, state, "get_set_str", gr_test_get_set_str, iters, test_flags);

    gr_test_iter(R, state, "one", gr_test_one, iters, test_flags);

    gr_test_iter(R, state, "mul: associative", gr_test_mul_associative, iters, test_flags);
/*
    if (gr_ctx_is_abelian_group(R) == T_TRUE)
        gr_test_iter(R, state, "mul: commutative", gr_test_mul_commutative, iters, test_flags);
*/
    gr_test_iter(R, state, "mul: aliasing", gr_test_mul_aliasing, iters, test_flags);

    gr_test_iter(R, state, "div: div then mul", gr_test_div_then_mul, iters, test_flags);
    gr_test_iter(R, state, "div: mul then div", gr_test_mul_then_div, iters, test_flags);

    gr_test_iter(R, state, "inv: multiplication", gr_test_inv_multiplication, iters, test_flags);
    gr_test_iter(R, state, "inv: involution", gr_test_inv_involution, iters, test_flags);

    gr_test_iter(R, state, "pow_ui: exponent addition", gr_test_pow_ui_exponent_addition, iters, test_flags);
    gr_test_iter(R, state, "pow_ui: aliasing", gr_test_pow_ui_exponent_addition, iters, test_flags);
    gr_test_iter(R, state, "pow_fmpz: exponent addition", gr_test_pow_fmpz_exponent_addition, iters, test_flags);

    gr_test_iter(R, state, "get_set_fexpr", gr_test_get_set_fexpr, iters, test_flags);

    flint_randclear(state);

    if (test_flags & GR_TEST_VERBOSE)
    {
        timeit_stop(timer);

        flint_printf("-------------------------------------------------------------------------------\n");
        flint_printf("Tests finished in %.3g cpu, %.3g wall\n", timer->cpu*0.001, timer->wall*0.001);
        flint_printf("===============================================================================\n\n");
    }
}
