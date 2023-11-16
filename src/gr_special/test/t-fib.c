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
test_fib_fmpz_rec1(flint_rand_t state)
{
    gr_ctx_t ctx;
    gr_ptr fx, fx1, fx2, fxfx1;
    fmpz_t x, x1, x2;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    GR_TMP_INIT4(fx, fx1, fx2, fxfx1, ctx);
    fmpz_init(x);
    fmpz_init(x1);
    fmpz_init(x2);

    if (gr_ctx_has_real_prec(ctx) == T_TRUE || gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
        fmpz_randtest(x, state, 100);
    else
        fmpz_randtest(x, state, 10);

    fmpz_add_ui(x1, x, 1);
    fmpz_add_ui(x2, x, 2);

    status |= gr_fib_fmpz(fx, x, ctx);
    status |= gr_fib_fmpz(fx1, x1, ctx);
    status |= gr_fib_fmpz(fx2, x2, ctx);
    status |= gr_add(fxfx1, fx, fx1, ctx);

    if (status == GR_SUCCESS && gr_equal(fxfx1, fx2, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        printf("x = "); fmpz_print(x); printf("\n");
        printf("x1 = "); fmpz_print(x1); printf("\n");
        printf("x2 = "); fmpz_print(x2); printf("\n");
        printf("fx = "); gr_println(fx, ctx); printf("\n");
        printf("fx1 = "); gr_println(fx1, ctx); printf("\n");
        printf("fxfx1 = "); gr_println(fxfx1, ctx); printf("\n");
        flint_abort();
    }

    GR_TMP_CLEAR4(fx, fx1, fx2, fxfx1, ctx);
    fmpz_clear(x);
    fmpz_clear(x1);
    fmpz_clear(x2);
    gr_ctx_clear(ctx);

    return status;
}

int
test_fib_ui_rec1(flint_rand_t state)
{
    gr_ctx_t ctx;
    gr_ptr fx, fx1, fx2, fxfx1;
    ulong x, x1, x2;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    GR_TMP_INIT4(fx, fx1, fx2, fxfx1, ctx);

    if (gr_ctx_has_real_prec(ctx) == T_TRUE || gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
    {
        x = n_randtest(state);
        x = FLINT_MIN(x, UWORD_MAX - 2);
    }
    else
    {
        x = n_randtest_bits(state, 10);
    }

    x1 = x + 1;
    x2 = x + 2;

    status |= gr_fib_ui(fx, x, ctx);
    status |= gr_fib_ui(fx1, x1, ctx);
    status |= gr_fib_ui(fx2, x2, ctx);
    status |= gr_add(fxfx1, fx, fx1, ctx);

    if (status == GR_SUCCESS && gr_equal(fxfx1, fx2, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        flint_printf("x = %wu", x); printf("\n");
        flint_printf("x1 = %wu", x1); printf("\n");
        flint_printf("x2 = %wu", x2); printf("\n");
        printf("fx = "); gr_println(fx, ctx); printf("\n");
        printf("fx1 = "); gr_println(fx1, ctx); printf("\n");
        printf("fxfx1 = "); gr_println(fxfx1, ctx); printf("\n");
        flint_abort();
    }

    GR_TMP_CLEAR4(fx, fx1, fx2, fxfx1, ctx);
    gr_ctx_clear(ctx);

    return status;
}

int
test_fib_vec(flint_rand_t state)
{
    gr_ctx_t ctx;
    slong len, i;
    gr_ptr x, xp;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    len = n_randint(state, 200);

    GR_TMP_INIT_VEC(xp, len, ctx);
    GR_TMP_INIT(x, ctx);

    status |= gr_fib_vec(xp, len, ctx);

    for (i = 0; i < len; i += 1 + n_randint(state, 10))
    {
        status |= gr_fib_ui(x, i, ctx);

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

TEST_FUNCTION_START(gr_special_fib, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
        test_fib_fmpz_rec1(state);

    for (iter = 0; iter < 1000; iter++)
        test_fib_ui_rec1(state);

    for (iter = 0; iter < 1000; iter++)
        test_fib_vec(state);

    TEST_FUNCTION_END(state);
}
