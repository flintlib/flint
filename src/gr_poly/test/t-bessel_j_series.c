/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_special.h"
#include "gr_poly.h"

/*
    Check consistency with
    besselj(nu, z) = z / 2nu * (besselj(nu - 1, z) + besselj(nu + 1, z)).
*/
int
test_bessel_j_series(flint_rand_t state)
{
    gr_ctx_t ctx;
    gr_ptr nu, tmp;
    gr_poly_t res1, res2, res3, z;
    slong len;
    int status = GR_SUCCESS;

    if (n_randint(state, 4) == 0)
	gr_ctx_init_real_arb(ctx, 128);
    else if (n_randint(state, 4) == 0)
	gr_ctx_init_complex_acb(ctx, 128);
    else
	gr_ctx_init_random_commutative_ring(ctx, state);

    GR_TMP_INIT2(nu, tmp, ctx);
    gr_poly_init(res1, ctx);
    gr_poly_init(res2, ctx);
    gr_poly_init(res3, ctx);
    gr_poly_init(z, ctx);

    len = n_randint(state, 10);

    status |= gr_randtest_not_zero(nu, state, ctx);
    GR_MUST_SUCCEED(gr_poly_randtest(z, state, n_randint(state, 20), ctx));

    if (n_randint(state, 4) != 0)
        status |= gr_poly_set_coeff_si(z, 0, 0, ctx);

    if (n_randint(state, 2))
    {
        status |= gr_poly_bessel_j_series(res1, nu, z, len, ctx);
    }
    else
    {
        /* test aliasing */
        status |= gr_poly_set(res1, z, ctx);
        status |= gr_poly_bessel_j_series(res1, nu, res1, len, ctx);
    }

    status |= gr_sub_si(tmp, nu, 1, ctx);
    status |= gr_poly_bessel_j_series(res2, tmp, z, len, ctx);

    status |= gr_add_si(tmp, nu, 1, ctx);
    status |= gr_poly_bessel_j_series(res3, tmp, z, len, ctx);

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_add(res2, res2, res3, ctx);
        status |= gr_poly_mullow(res2, res2, z, len, ctx);
        status |= gr_mul_si(tmp, nu, 2, ctx);
        status |= gr_poly_div_scalar(res2, res2, tmp, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(res1, res2, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("len = %{slong}\n", len);
            printf("nu = "); gr_println(nu, ctx); printf("\n\n");
            printf("z = "); gr_poly_print(z, ctx); printf("\n\n");
            flint_abort();
        }
    }

    GR_TMP_CLEAR2(nu, tmp, ctx);
    gr_poly_clear(res1, ctx);
    gr_poly_clear(res2, ctx);
    gr_poly_clear(res3, ctx);
    gr_poly_clear(z, ctx);

    gr_ctx_clear(ctx);

    return status;
}

/*
    Check that zero input gives Bessel function at zero.
*/
int
test_bessel_j_series_zero(flint_rand_t state)
{
    gr_ctx_t ctx;
    gr_ptr nu, besselj, coeff;
    gr_poly_t res, z;
    slong len;
    int status = GR_SUCCESS;

    if (n_randint(state, 4) == 0)
	gr_ctx_init_real_arb(ctx, 128);
    else if (n_randint(state, 4) == 0)
	gr_ctx_init_complex_acb(ctx, 128);
    else
	gr_ctx_init_random_commutative_ring(ctx, state);

    GR_TMP_INIT3(nu, besselj, coeff, ctx);
    gr_poly_init(res, ctx);
    gr_poly_init(z, ctx);

    len = n_randint(state, 10);

    status |= gr_randtest(nu, state, ctx);

    status |= gr_poly_bessel_j_series(res, nu, z, len, ctx);
    status |= gr_bessel_j(besselj, nu, besselj, ctx);

    if (status == GR_SUCCESS)
    {
        if (gr_poly_is_scalar(res, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("len = %{slong}\n", len);
            printf("nu = "); gr_println(nu, ctx); printf("\n\n");
            printf("z = "); gr_poly_print(z, ctx); printf("\n\n");
            printf("res = "); gr_poly_print(res, ctx); printf("\n\n");
            flint_abort();
        }

        status |= gr_poly_get_coeff_scalar(coeff, res, 0, ctx);

        if (len > 0 && gr_equal(coeff, besselj, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("len = %{slong}\n", len);
            printf("nu = "); gr_println(nu, ctx); printf("\n\n");
            printf("z = "); gr_poly_print(z, ctx); printf("\n\n");
            printf("coeff = "); gr_println(coeff, ctx); printf("\n\n");
            printf("besselj = "); gr_println(besselj, ctx); printf("\n\n");
            flint_abort();
        }
    }

    GR_TMP_CLEAR3(nu, besselj, coeff, ctx);
    gr_poly_clear(res, ctx);
    gr_poly_clear(z, ctx);

    gr_ctx_clear(ctx);

    return status;
}


TEST_FUNCTION_START(gr_poly_bessel_j_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        //test_bessel_j_series(state);
    }

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        test_bessel_j_series_zero(state);
    }


    TEST_FUNCTION_END(state);
}
