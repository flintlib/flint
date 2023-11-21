/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr_special.h"

int
test_fac_rec1(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    gr_ptr x, x1, fx, fx1, x1fx;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    GR_TMP_INIT5(x, x1, fx, fx1, x1fx, ctx);

    if (gr_ctx_is_exact(ctx) == T_TRUE)
        status |= gr_set_si(x, n_randint(state, 200) - 100, ctx);
    else
        status |= gr_randtest(x, state, ctx);

    status |= gr_add_ui(x1, x, 1, ctx);

    if (which == 0)
    {
        status |= gr_fac(fx, x, ctx);
        status |= gr_fac(fx1, x1, ctx);
        status |= gr_mul(x1fx, fx, x1, ctx);
    }
    else
    {
        status |= gr_rfac(fx, x, ctx);
        status |= gr_rfac(fx1, x1, ctx);
        status |= gr_div(x1fx, fx, x1, ctx);
    }

    if (status == GR_SUCCESS && gr_equal(fx1, x1fx, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        printf("x = "); gr_println(x, ctx); printf("\n");
        printf("x1 = "); gr_println(x1, ctx); printf("\n");
        printf("fx = "); gr_println(fx, ctx); printf("\n");
        printf("fx1 = "); gr_println(fx1, ctx); printf("\n");
        printf("x1fx = "); gr_println(x1fx, ctx); printf("\n");
        flint_abort();
    }

    GR_TMP_CLEAR5(x, x1, fx, fx1, x1fx, ctx);
    gr_ctx_clear(ctx);

    return status;
}

int
test_fac_fmpz_rec1(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    gr_ptr fx, fx1, x1fx;
    fmpz_t x, x1;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    GR_TMP_INIT3(fx, fx1, x1fx, ctx);
    fmpz_init(x);
    fmpz_init(x1);

    if (gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        fmpz_randtest(x, state, 100);
    }
    else
    {
        if (gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
            fmpz_randtest(x, state, 12);
        else
            fmpz_randtest(x, state, 8);
    }

    fmpz_add_ui(x1, x, 1);

    if (which == 0)
    {
        status |= gr_fac_fmpz(fx, x, ctx);
        status |= gr_fac_fmpz(fx1, x1, ctx);
        status |= gr_mul_fmpz(x1fx, fx, x1, ctx);
    }
    else
    {
        status |= gr_rfac_fmpz(fx, x, ctx);
        status |= gr_rfac_fmpz(fx1, x1, ctx);
        status |= gr_div_fmpz(x1fx, fx, x1, ctx);
    }

    if (status == GR_SUCCESS && gr_equal(fx1, x1fx, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        printf("x = "); fmpz_print(x); printf("\n");
        printf("x1 = "); fmpz_print(x1); printf("\n");
        printf("fx = "); gr_println(fx, ctx); printf("\n");
        printf("fx1 = "); gr_println(fx1, ctx); printf("\n");
        printf("x1fx = "); gr_println(x1fx, ctx); printf("\n");
        flint_abort();
    }

    GR_TMP_CLEAR3(fx, fx1, x1fx, ctx);
    fmpz_clear(x);
    fmpz_clear(x1);
    gr_ctx_clear(ctx);

    return status;
}

int
test_fac_ui_rec1(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    gr_ptr fx, fx1, x1fx;
    ulong x, x1;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    GR_TMP_INIT3(fx, fx1, x1fx, ctx);

    if (gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        x = n_randtest(state);
        x -= (x == UWORD_MAX);
    }
    else
    {
        if (gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
            x = n_randtest_bits(state, 12);
        else
            x = n_randtest_bits(state, 8);
    }

    x1 = x + 1;

    if (which == 0)
    {
        status |= gr_fac_ui(fx, x, ctx);
        status |= gr_fac_ui(fx1, x1, ctx);
        status |= gr_mul_ui(x1fx, fx, x1, ctx);
    }
    else
    {
        status |= gr_rfac_ui(fx, x, ctx);
        status |= gr_rfac_ui(fx1, x1, ctx);
        status |= gr_div_ui(x1fx, fx, x1, ctx);
    }

    if (status == GR_SUCCESS && gr_equal(fx1, x1fx, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        flint_printf("x = %wu", x); printf("\n");
        flint_printf("x1 = %wu", x1); printf("\n");
        printf("fx = "); gr_println(fx, ctx); printf("\n");
        printf("fx1 = "); gr_println(fx1, ctx); printf("\n");
        printf("x1fx = "); gr_println(x1fx, ctx); printf("\n");
        flint_abort();
    }

    GR_TMP_CLEAR3(fx, fx1, x1fx, ctx);
    gr_ctx_clear(ctx);

    return status;
}

int
test_fac_vec(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    slong len, i;
    gr_ptr x, xp;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    len = n_randint(state, 100);

    GR_TMP_INIT_VEC(xp, len, ctx);
    GR_TMP_INIT(x, ctx);

    if (which == 0)
        status |= gr_fac_vec(xp, len, ctx);
    else
        status |= gr_rfac_vec(xp, len, ctx);

    for (i = 0; i < len; i += 1 + n_randint(state, 10))
    {
        if (which == 0)
            status |= gr_fac_ui(x, i, ctx);
        else
            status |= gr_rfac_ui(x, i, ctx);

        if (status == GR_SUCCESS && gr_equal(x, GR_ENTRY(xp, i, ctx->sizeof_elem), ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            printf("x = "); gr_println(x, ctx); printf("\n");
            printf("xp = "); gr_println(GR_ENTRY(xp, i, ctx->sizeof_elem), ctx); printf("\n");
            flint_abort();
        }
    }

    GR_TMP_CLEAR(x, ctx);
    GR_TMP_CLEAR_VEC(xp, len, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_special_fac, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test_fac_rec1(state, 0);
        test_fac_rec1(state, 1);
    }

    for (iter = 0; iter < 1000; iter++)
    {
        test_fac_fmpz_rec1(state, 0);
        test_fac_fmpz_rec1(state, 1);
    }

    for (iter = 0; iter < 1000; iter++)
    {
        test_fac_ui_rec1(state, 0);
        test_fac_ui_rec1(state, 1);
    }

    for (iter = 0; iter < 1000; iter++)
    {
        test_fac_vec(state, 0);
        test_fac_vec(state, 1);
    }

    TEST_FUNCTION_END(state);
}
