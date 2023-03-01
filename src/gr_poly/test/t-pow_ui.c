/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
test(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    /* Test A^n * A^m = A^(n+m) */
    gr_poly_t A, An, Am, AnAm, Anm;
    ulong n, m;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(A, ctx);
    gr_poly_init(An, ctx);
    gr_poly_init(Am, ctx);
    gr_poly_init(AnAm, ctx);
    gr_poly_init(Anm, ctx);

    n = n_randint(state, 5);
    m = n_randint(state, 5);

    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 8, ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(An, state, 8, ctx));

    if (which == 0)
    {
        status |= gr_poly_pow_ui(An, A, n, ctx);
        status |= gr_poly_set(Am, A, ctx);
        status |= gr_poly_pow_ui(Am, Am, m, ctx);
        status |= gr_poly_mul(AnAm, An, Am, ctx);
        status |= gr_poly_pow_ui(Anm, A, n + m, ctx);
    }
    else
    {
        status |= gr_poly_pow_ui_binexp(An, A, n, ctx);
        status |= gr_poly_set(Am, A, ctx);
        status |= gr_poly_pow_ui_binexp(Am, Am, m, ctx);
        status |= gr_poly_mul(AnAm, An, Am, ctx);
        status |= gr_poly_pow_ui_binexp(Anm, A, n + m, ctx);
    }

    if (gr_poly_equal(AnAm, Anm, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n\n");
        gr_ctx_println(ctx);
        flint_printf("n = %wu, m = %wu\n", n, m);
        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
        flint_printf("An = "); gr_poly_print(An, ctx); flint_printf("\n");
        flint_printf("Am = "); gr_poly_print(Am, ctx); flint_printf("\n");
        flint_printf("AnAm = "); gr_poly_print(AnAm, ctx); flint_printf("\n");
        flint_printf("Anm = "); gr_poly_print(Anm, ctx); flint_printf("\n");
        flint_abort();
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(An, ctx);
    gr_poly_clear(Am, ctx);
    gr_poly_clear(AnAm, ctx);
    gr_poly_clear(Anm, ctx);

    gr_ctx_clear(ctx);

    return status;
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("pow_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        test(state, n_randint(state, 2));
    }

    flint_randclear(state);
    flint_cleanup_master();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
