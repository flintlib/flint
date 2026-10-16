/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "gr_poly.h"
#include "long_extras.h"

/*
    Check that the Taylor coefficients satisfy the recurrence from
    https://fungrim.org/entry/9b2f38/
*/
int
test_recurrence(flint_rand_t state)
{
    gr_ctx_t ctx;
    gr_ptr res, nu, z, x, tmp;
    slong sz;
    int status = GR_SUCCESS;
    int len;

    if (n_randint(state, 4) == 0)
	gr_ctx_init_real_arb(ctx, 128);
    else if (n_randint(state, 4) == 0)
	gr_ctx_init_complex_acb(ctx, 128);
    else
	gr_ctx_init_random_commutative_ring(ctx, state);

    sz = ctx->sizeof_elem;

    GR_TMP_INIT4(nu, z, x, tmp, ctx);

    len = n_randint(state, 10);
    res = gr_heap_init_vec(len, ctx);

    if (n_randint(state, 4) == 0)
	status |= gr_set_si(nu, z_randint(state, len + 1), ctx);
    else
	GR_MUST_SUCCEED(gr_randtest(nu, state, ctx));

    if (n_randint(state, 4) == 0)
	GR_MUST_SUCCEED(gr_zero(z, ctx));
    else
	GR_MUST_SUCCEED(gr_randtest(z, state, ctx));

    status |= gr_bessel_j_jet(res, nu, z, len, ctx);

    if (status == GR_SUCCESS)
    {
	slong r;

	if (len > 0)
	{
	    status |= gr_bessel_j(tmp, nu, z, ctx);

	    if (status == GR_SUCCESS && gr_equal(GR_ENTRY(res, 0, sz), tmp, ctx) == T_FALSE)
	    {
		flint_printf("FAIL\n");
		printf("ctx = "); gr_ctx_println(ctx); printf("\n");
		printf("nu = "); gr_println(nu, ctx); printf("\n");
		printf("z = "); gr_println(z, ctx); printf("\n");
		flint_printf("len = %{slong}\n", len);
		flint_abort();
	    }
	}

	for (r = 0; r < len - 4; r++)
	{
	    status |= gr_set(x, GR_ENTRY(res, r, sz), ctx);

	    status |= gr_mul(tmp, GR_ENTRY(res, r + 1, sz), z, ctx);
	    status |= gr_mul_si(tmp, tmp, 2, ctx);
	    status |= gr_add(x, x, tmp, ctx);

	    status |= gr_sqr(tmp, z, ctx);
	    status |= gr_submul(tmp, nu, nu, ctx);
	    status |= gr_add_si(tmp, tmp, r * (r + 4) + 4, ctx);
	    status |= gr_mul(tmp, tmp, GR_ENTRY(res, r + 2, sz), ctx);
	    status |= gr_add(x, x, tmp, ctx);

	    status |= gr_mul(tmp, GR_ENTRY(res, r + 3, sz), z, ctx);
	    status |= gr_mul_si(tmp, tmp, 2 * r * r + 11 * r + 15, ctx);
	    status |= gr_add(x, x, tmp, ctx);

	    status |= gr_sqr(tmp, z, ctx);
	    status |= gr_mul(tmp, tmp, GR_ENTRY(res, r + 4, sz), ctx);
	    status |= gr_mul_si(tmp, tmp, r * r + 7 * r + 12, ctx);
	    status |= gr_add(x, x, tmp, ctx);

	    if (gr_is_zero(x, ctx) == T_FALSE)
	    {
		flint_printf("FAIL\n");
		printf("ctx = "); gr_ctx_println(ctx); printf("\n");
		printf("nu = "); gr_println(nu, ctx); printf("\n");
		printf("z = "); gr_println(z, ctx); printf("\n");
		flint_printf("len = %{slong}\n", len);
		flint_printf("r = %{slong}\n", r);
		printf("x = "); gr_println(x, ctx); printf("\n");
		flint_abort();
	    }
	}
    }

    GR_TMP_CLEAR4(nu, z, x, tmp, ctx);
    gr_heap_clear_vec(res, len, ctx);
    gr_ctx_clear(ctx);

    return status;
}

/*
    Check consistency with
    besselj(nu, z) = z / 2nu * (besselj(nu - 1, z) + besselj(nu + 1, z)).
*/
int
test_consistency_with_bessel_j(flint_rand_t state)
{
    gr_ctx_t ctx;
    gr_ptr res1, res2, res3, nu, z, tmp;
    slong sz;
    int status = GR_SUCCESS;
    int len;

    if (n_randint(state, 4) == 0)
	gr_ctx_init_real_arb(ctx, 128);
    else if (n_randint(state, 4) == 0)
	gr_ctx_init_complex_acb(ctx, 128);
    else
	gr_ctx_init_random_commutative_ring(ctx, state);

    sz = ctx->sizeof_elem;

    GR_TMP_INIT2(nu, z, ctx);

    len = 1 + n_randint(state, 10);
    res1 = gr_heap_init_vec(len, ctx);
    res2 = gr_heap_init_vec(len, ctx);
    res3 = gr_heap_init_vec(len, ctx);
    tmp = gr_heap_init_vec(len, ctx);

    if (n_randint(state, 4) == 0)
	status |= gr_set_si(nu, z_randint(state, len + 1), ctx);
    else
	GR_MUST_SUCCEED(gr_randtest(nu, state, ctx));

    if (n_randint(state, 4) == 0)
	GR_MUST_SUCCEED(gr_zero(z, ctx));
    else
	GR_MUST_SUCCEED(gr_randtest(z, state, ctx));

    status |= gr_bessel_j_jet(res1, nu, z, len, ctx);

    status |= gr_sub_si(tmp, nu, 1, ctx);
    status |= gr_bessel_j_jet(res2, tmp, z, len, ctx);

    status |= gr_add_si(tmp, nu, 1, ctx);
    status |= gr_bessel_j_jet(res3, tmp, z, len, ctx);

    if (status == GR_SUCCESS)
    {
	status |= _gr_vec_add(res2, res2, res3, len, ctx);

	status |= gr_set(GR_ENTRY(tmp, 0, sz), z, ctx);
	if (len > 1)
	    status |= gr_one(GR_ENTRY(tmp, 1, sz), ctx);
	status |= _gr_vec_div_scalar(tmp, tmp, len, nu, ctx);
	status |= _gr_vec_mul_scalar_2exp_si(tmp, tmp, len, -1, ctx);

	status |= _gr_poly_mullow(res3, res2, len, tmp, len, len, ctx);

	if (status == GR_SUCCESS && _gr_vec_equal(res1, res3, len, ctx) == T_FALSE)
	{
	    flint_printf("FAIL\n");
	    printf("ctx = "); gr_ctx_println(ctx); printf("\n");
	    printf("nu = "); gr_println(nu, ctx); printf("\n");
	    printf("z = "); gr_println(z, ctx); printf("\n");
	    flint_printf("len = %{slong}\n", len);
	    flint_abort();
	}
    }

    GR_TMP_CLEAR2(nu, z, ctx);
    gr_heap_clear_vec(res1, len, ctx);
    gr_heap_clear_vec(res2, len, ctx);
    gr_heap_clear_vec(res3, len, ctx);
    gr_heap_clear_vec(tmp, len, ctx);
    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_bessel_j_jet, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
        test_recurrence(state);

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
        test_consistency_with_bessel_j(state);

    TEST_FUNCTION_END(state);
}
