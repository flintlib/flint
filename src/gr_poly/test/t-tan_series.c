/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_special.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

static int
gr_poly_tan_series_ref(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status;
    gr_poly_t s, c;
    gr_poly_init(s, ctx);
    gr_poly_init(c, ctx);
    status = gr_poly_sin_cos_series(s, c, h, n, ctx);
    status |= gr_poly_div_series(f, s, c, n, ctx);
    gr_poly_clear(s, ctx);
    gr_poly_clear(c, ctx);
    return status;
}

static int
_gr_poly_tanh_series_ref(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    flen = FLINT_MIN(flen, n);

    gr_ptr t, u;
    GR_TMP_INIT_VEC(t, n, ctx);
    GR_TMP_INIT_VEC(u, n, ctx);

    status |= _gr_vec_mul_scalar_2exp_si(t, f, flen, 1, ctx);
    status |= _gr_poly_exp_series(u, t, flen, n, ctx);
    status |= _gr_vec_set(t, u, n, ctx);
    status |= gr_sub_ui(t, t, 1, ctx);
    status |= gr_add_ui(u, u, 1, ctx);
    status |= _gr_poly_div_series(res, t, n, u, n, n, ctx);

    GR_TMP_CLEAR_VEC(t, n, ctx);
    GR_TMP_CLEAR_VEC(u, n, ctx);

    return status;
}

static int
gr_poly_tanh_series_ref(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen = h->length;

    if (n == 0 || hlen == 0)
        return gr_poly_zero(f, ctx);

    gr_poly_fit_length(f, n, ctx);
    status |= _gr_poly_tanh_series_ref(f->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(f, n, ctx);
    _gr_poly_normalise(f, ctx);
    return status;
}

static int
gr_poly_cot_series_ref(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status;
    status = gr_poly_tan_series_ref(f, h, n, ctx);
    status |= gr_poly_inv_series(f, f, n, ctx);
    return status;
}

static int
gr_poly_coth_series_ref(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status;
    status = gr_poly_tanh_series_ref(f, h, n, ctx);
    status |= gr_poly_inv_series(f, f, n, ctx);
    return status;
}

static int
gr_poly_tan_pi_series_ref(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t = gr_heap_init(ctx);
    status |= gr_pi(t, ctx);
    status |= gr_poly_mul_scalar(f, h, t, ctx);
    status |= gr_poly_tan_series_ref(f, f, n, ctx);
    gr_heap_clear(t, ctx);
    return status;
}

static int
gr_poly_cot_pi_series_ref(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t = gr_heap_init(ctx);
    status |= gr_pi(t, ctx);
    status |= gr_poly_mul_scalar(f, h, t, ctx);
    status |= gr_poly_cot_series_ref(f, f, n, ctx);
    gr_heap_clear(t, ctx);
    return status;
}

TEST_FUNCTION_START(gr_poly_tan_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong n;
        gr_poly_t a, T, Tref;
        int status = GR_SUCCESS;
        int func, aliasing, algorithm;

        if (n_randint(state, 4) == 0)
            gr_ctx_init_real_arb(ctx, 128);
        else if (n_randint(state, 4) == 0)
            gr_ctx_init_complex_acb(ctx, 128);
        else
            gr_ctx_init_random_commutative_ring(ctx, state);

        gr_poly_init(a, ctx);
        gr_poly_init(T, ctx);
        gr_poly_init(Tref, ctx);

        if (ctx->methods == _ca_methods)
            n = 1 + n_randint(state, 4);
        else if (gr_ctx_is_finite(ctx) == T_TRUE)
            n = 1 + n_randint(state, 80);
        else
            n = 1 + n_randint(state, 30);

        GR_MUST_SUCCEED(gr_poly_randtest(a, state, n + 2, ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(Tref, state, n + 2, ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(T, state, n + 2, ctx));

        if (n_randint(state, 4) != 0)
            status |= gr_poly_set_coeff_si(a, 0, 0, ctx);

        do {
            func = n_randint(state, 8);
        } while (func == 5 || func == 7);

        /* Test basecase vs algorithm */

        status |= gr_poly_tan_series_basecase(Tref, a, n + n_randint(state, 2), func, ctx);
        status |= gr_poly_truncate(Tref, Tref, n, ctx);

        aliasing = n_randint(state, 2);

        algorithm = n_randint(state, 3);

        if (algorithm == 0)
        {
            if (aliasing)
            {
                status |= gr_poly_set(T, a, ctx);
                status |= gr_poly_tan_series_newton(T, T, n, n_randint(state, 20), func, ctx);
            }
            else
            {
                status |= gr_poly_tan_series_newton(T, a, n, n_randint(state, 20), func, ctx);
            }
        }
        else if (algorithm == 1)
        {
            if (aliasing)
            {
                status |= gr_poly_set(T, a, ctx);
                status |= gr_poly_tan_series_sine_cosine(T, T, n, func, ctx);
            }
            else
            {
                status |= gr_poly_tan_series_sine_cosine(T, a, n, func, ctx);
            }
        }
        else if (algorithm == 2)
        {
            if (aliasing)
            {
                status |= gr_poly_set(T, a, ctx);
                status |= gr_poly_tan_series_exponential(T, T, n, func, ctx);
            }
            else
            {
                status |= gr_poly_tan_series_exponential(T, a, n, func, ctx);
            }
        }

        if (status == GR_SUCCESS && (gr_poly_equal(T, Tref, ctx) == T_FALSE))
        {
            flint_printf("FAIL: gr_poly_tan_series (basecase vs algorithm %d)\n\n", algorithm);
            gr_ctx_println(ctx);
            flint_printf("func = %d, aliasing = %d\n", func, aliasing);
            flint_printf("n = %wd\n", n);
            flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("Tref = "); gr_poly_print(Tref, ctx); flint_printf("\n\n");
            flint_printf("T    = "); gr_poly_print(T, ctx); flint_printf("\n\n");
            flint_abort();
        }

        /* Test default function vs ref */
        status = 0;

        switch (func)
        {
            case 0: status |= gr_poly_tan_series_ref(Tref, a, n, ctx); break;
            case 1: status |= gr_poly_tanh_series_ref(Tref, a, n, ctx); break;
            case 2: status |= gr_poly_cot_series_ref(Tref, a, n, ctx); break;
            case 3: status |= gr_poly_coth_series_ref(Tref, a, n, ctx); break;
            case 4: status |= gr_poly_tan_pi_series_ref(Tref, a, n, ctx); break;
            case 6: status |= gr_poly_cot_pi_series_ref(Tref, a, n, ctx); break;
            default: flint_abort();
        }

        if (aliasing)
        {
            status |= gr_poly_set(T, a, ctx);

            switch (func)
            {
                case 0: status |= gr_poly_tan_series(T, T, n, ctx); break;
                case 1: status |= gr_poly_tanh_series(T, T, n, ctx); break;
                case 2: status |= gr_poly_cot_series(T, T, n, ctx); break;
                case 3: status |= gr_poly_coth_series(T, T, n, ctx); break;
                case 4: status |= gr_poly_tan_pi_series(T, T, n, ctx); break;
                case 6: status |= gr_poly_cot_pi_series(T, T, n, ctx); break;
                default: flint_abort();
            }
        }
        else
        {
            switch (func)
            {
                case 0: status |= gr_poly_tan_series(T, a, n, ctx); break;
                case 1: status |= gr_poly_tanh_series(T, a, n, ctx); break;
                case 2: status |= gr_poly_cot_series(T, a, n, ctx); break;
                case 3: status |= gr_poly_coth_series(T, a, n, ctx); break;
                case 4: status |= gr_poly_tan_pi_series(T, a, n, ctx); break;
                case 6: status |= gr_poly_cot_pi_series(T, a, n, ctx); break;
                default: flint_abort();
            }
        }

        if (status == GR_SUCCESS && (gr_poly_equal(T, Tref, ctx) == T_FALSE))
        {
            flint_printf("FAIL: gr_poly_tan_series (ref vs default)\n\n");
            gr_ctx_println(ctx);
            flint_printf("func = %d, aliasing = %d\n", func, aliasing);
            flint_printf("n = %wd\n", n);
            flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("Tref = "); gr_poly_print(Tref, ctx); flint_printf("\n\n");
            flint_printf("T    = "); gr_poly_print(T, ctx); flint_printf("\n\n");
            flint_abort();
        }

        gr_poly_clear(a, ctx);
        gr_poly_clear(Tref, ctx);
        gr_poly_clear(T, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

