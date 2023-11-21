/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "gr_special.h"

int
test_chebyshev_fmpz_rec1(flint_rand_t state)
{
    gr_ctx_t ctx;
    gr_ptr ft, ft1, ft2, t, g;
    fmpz_t n, n1, n2;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    GR_TMP_INIT5(ft, ft1, ft2, g, t, ctx);
    fmpz_init(n);
    fmpz_init(n1);
    fmpz_init(n2);

    if (gr_ctx_has_real_prec(ctx) == T_TRUE || gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
        fmpz_randtest(n, state, 100);
    else
        fmpz_randtest(n, state, 5);

    fmpz_add_ui(n1, n, 1);
    fmpz_add_ui(n2, n, 2);

    status |= gr_randtest(t, state, ctx);

    if (n_randint(state, 2))
    {
        status |= gr_generic_chebyshev_t_fmpz(ft, n, t, ctx);
        status |= gr_generic_chebyshev_t_fmpz(ft1, n1, t, ctx);
        status |= gr_generic_chebyshev_t_fmpz(ft2, n2, t, ctx);
    }
    else
    {
        status |= gr_generic_chebyshev_u_fmpz(ft, n, t, ctx);
        status |= gr_generic_chebyshev_u_fmpz(ft1, n1, t, ctx);
        status |= gr_generic_chebyshev_u_fmpz(ft2, n2, t, ctx);
    }

    status |= gr_mul(g, ft1, t, ctx);
    status |= gr_mul_two(g, g, ctx);
    status |= gr_sub(g, g, ft, ctx);

    if (status == GR_SUCCESS && gr_equal(g, ft2, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        printf("n = "); fmpz_print(n); printf("\n");
        printf("n1 = "); fmpz_print(n1); printf("\n");
        printf("n2 = "); fmpz_print(n2); printf("\n");
        printf("t = "); gr_println(t, ctx); printf("\n");
        printf("ft = "); gr_println(ft, ctx); printf("\n");
        printf("ft1 = "); gr_println(ft1, ctx); printf("\n");
        printf("ft2 = "); gr_println(ft2, ctx); printf("\n");
        printf("g = "); gr_println(g, ctx); printf("\n");
        flint_abort();
    }

    GR_TMP_CLEAR5(ft, ft1, ft2, g, t, ctx);
    fmpz_clear(n);
    fmpz_clear(n1);
    fmpz_clear(n2);
    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_special_chebyshev, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
        test_chebyshev_fmpz_rec1(state);

    TEST_FUNCTION_END(state);
}
