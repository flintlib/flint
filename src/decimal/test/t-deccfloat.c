/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "decimal.h"
#include "gr.h"
#include "gr_poly.h"

/* parse s in ctx and check that the result prints as expected */
static void
_ccheck_parse(const char * s, const char * expected, gr_ctx_t ctx)
{
    deccfloat_t x;
    char * t;
    int status;

    deccfloat_init(x, ctx);
    status = gr_set_str(x, s, ctx);
    t = (status == GR_SUCCESS) ? deccfloat_get_str(x, ctx) : NULL;

    if (status != GR_SUCCESS || strcmp(t, expected) != 0)
    {
        flint_printf("FAIL: parse\n");
        gr_ctx_println(ctx);
        flint_printf("s = %s\n", s);
        flint_printf("expected = %s\n", expected);
        flint_printf("got = %s (status %d)\n", t ? t : "(null)", status);
        flint_abort();
    }

    flint_free(t);
    deccfloat_clear(x, ctx);
}

TEST_FUNCTION_START(deccfloat, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;

        gr_ctx_init_deccfloat_randtest(ctx, state, 60);

        if (gr_ctx_is_ring(ctx) == T_TRUE)
            gr_test_ring(ctx, 50, 0);
        else
            gr_test_floating_point(ctx, 50, 0);

        gr_ctx_clear(ctx);
    }

    {
        gr_ctx_t ctx;
        slong prec;

        for (prec = 1; prec <= 100; prec = 2 * prec + 1)
        {
            gr_ctx_init_deccfloat(ctx, prec, DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN);
            gr_test_floating_point(ctx, 20 * flint_test_multiplier(), 0);
            gr_ctx_clear(ctx);

            gr_ctx_init_deccfloat(ctx, prec, 0);
            decimal_ctx_set_rnd_im(ctx, DECIMAL_RND_FLOOR);
            gr_test_floating_point(ctx, 20 * flint_test_multiplier(), 0);
            gr_ctx_clear(ctx);
        }

        gr_ctx_init_deccfloat(ctx, DECIMAL_PREC_EXACT, 0);
        gr_test_ring(ctx, 100 * flint_test_multiplier(), 0);
        gr_ctx_clear(ctx);
    }

    /* strings */
    {
        gr_ctx_t ctx;

        gr_ctx_init_deccfloat(ctx, 20, DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN);

        _ccheck_parse("1.5", "1.5", ctx);
        _ccheck_parse("I", "1*I", ctx);
        _ccheck_parse("-I", "-1*I", ctx);
        _ccheck_parse("2*I", "2*I", ctx);
        _ccheck_parse("1+2*I", "(1 + 2*I)", ctx);
        _ccheck_parse("(1 + 2*I)", "(1 + 2*I)", ctx);
        _ccheck_parse("(1 - 2*I)", "(1 - 2*I)", ctx);
        _ccheck_parse("(-1.5e-7 - 2.25e10*I)", "(-1.5e-7 - 22500000000*I)", ctx);
        _ccheck_parse("(1 + 2*I) * (3 - 4*I)", "(11 + 2*I)", ctx);
        _ccheck_parse("(1 + 2*I) / (3 - 4*I)", "(-0.2 + 0.4*I)", ctx);
        _ccheck_parse("I^2", "-1", ctx);
        _ccheck_parse("sqrt(-4)", "2*I", ctx);
        _ccheck_parse("(3+4*I)^(1/2)", "(2 + 1*I)", ctx);
        _ccheck_parse("1/3 + I/3", "(0.33333333333333333333 + 0.33333333333333333333*I)", ctx);
        _ccheck_parse("inf", "inf", ctx);
        _ccheck_parse("-inf*I", "-inf*I", ctx);
        _ccheck_parse("(1 + nan*I)", "(nan + nan*I)", ctx);

        gr_ctx_clear(ctx);
    }

    /* string roundtrips, also in polynomials */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, pctx;
        deccfloat_t x, y;
        gr_poly_t f, g;
        char * s;
        int status;

        gr_ctx_init_deccfloat_randtest(ctx, state, 30);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        gr_ctx_init_gr_poly(pctx, ctx);

        deccfloat_init(x, ctx);
        deccfloat_init(y, ctx);

        GR_MUST_SUCCEED(deccfloat_randtest(x, state, ctx));
        s = deccfloat_get_str(x, ctx);
        status = gr_set_str(y, s, ctx);

        if (status != GR_SUCCESS || (deccfloat_equal(x, y, ctx) != T_TRUE && !_deccfloat_is_nan(x)))
        {
            flint_printf("FAIL: string roundtrip\n");
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}\n", x, ctx);
            flint_printf("s = %s\n", s);
            flint_printf("y = %{gr} (%d)\n", y, ctx, status);
            flint_abort();
        }

        flint_free(s);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        status = gr_poly_randtest(f, state, 1 + n_randint(state, 5), ctx);
        /* polynomial arithmetic is not meaningful with infinite coefficients */
        {
            slong i;
            for (i = 0; i < f->length; i++)
                if (!_deccfloat_is_finite(gr_poly_coeff_ptr(f, i, ctx)))
                    GR_MUST_SUCCEED(deccfloat_zero(gr_poly_coeff_ptr(f, i, ctx), ctx));
            _gr_poly_normalise(f, ctx);
        }
        status |= gr_get_str(&s, f, pctx);
        status |= gr_set_str(g, s, pctx);

        if (status != GR_SUCCESS || (gr_poly_equal(f, g, ctx) != T_TRUE && gr_poly_equal(f, g, ctx) != T_UNKNOWN))
        {
            flint_printf("FAIL: polynomial roundtrip\n");
            gr_ctx_println(ctx);
            flint_printf("f = %s\n", s);
            flint_printf("g = %{gr}  (%d)\n", g, pctx, status);
            flint_abort();
        }

        flint_free(s);
        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        deccfloat_clear(x, ctx);
        deccfloat_clear(y, ctx);
        gr_ctx_clear(pctx);
        gr_ctx_clear(ctx);
    }

    /* conversions between decimal contexts */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx1, ctx2, rctx;
        deccfloat_t x, y, z;
        decfloat_t r;
        int status;

        gr_ctx_init_deccfloat_randtest(ctx1, state, 30);
        gr_ctx_init_deccfloat_randtest(ctx2, state, 30);
        gr_ctx_init_decfloat_randtest(rctx, state, 30);
        decimal_ctx_set_exp_limits(ctx1, WORD_MIN, WORD_MAX);
        decimal_ctx_set_exp_limits(ctx2, WORD_MIN, WORD_MAX);
        decimal_ctx_set_exp_limits(rctx, WORD_MIN, WORD_MAX);
        decimal_ctx_set_prec(ctx2, DECIMAL_PREC_EXACT);
        decimal_ctx_set_prec(rctx, DECIMAL_PREC_EXACT);
        decimal_ctx_set_flags(ctx2, DECIMAL_CTX_FLAGS(ctx2) | DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN);
        decimal_ctx_set_flags(rctx, DECIMAL_CTX_FLAGS(rctx) | DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN);

        deccfloat_init(x, ctx1);
        deccfloat_init(y, ctx2);
        deccfloat_init(z, ctx1);
        decfloat_init(r, rctx);

        GR_MUST_SUCCEED(deccfloat_randtest(x, state, ctx1));

        /* exact conversion to a context with exact precision and back */
        status = gr_set_other(y, x, ctx1, ctx2);
        status |= gr_set_other(z, y, ctx2, ctx1);

        if (status != GR_SUCCESS || (deccfloat_equal(x, z, ctx1) != T_TRUE && !_deccfloat_is_nan(x)))
        {
            flint_printf("FAIL: context roundtrip\n");
            gr_ctx_println(ctx1);
            gr_ctx_println(ctx2);
            flint_printf("x = %{gr}\n", x, ctx1);
            flint_printf("y = %{gr}\n", y, ctx2);
            flint_printf("z = %{gr} (%d)\n", z, ctx1, status);
            flint_abort();
        }

        /* real numbers convert to and from the real type */
        decfloat_zero(&x->im, ctx1);
        status = gr_set_other(r, x, ctx1, rctx);
        status |= gr_set_other(z, r, rctx, ctx1);

        if (status != GR_SUCCESS || (deccfloat_equal(x, z, ctx1) != T_TRUE && !_deccfloat_is_nan(x)))
        {
            flint_printf("FAIL: real roundtrip\n");
            gr_ctx_println(ctx1);
            gr_ctx_println(rctx);
            flint_printf("x = %{gr}\n", x, ctx1);
            flint_printf("r = %{gr}\n", r, rctx);
            flint_printf("z = %{gr} (%d)\n", z, ctx1, status);
            flint_abort();
        }

        deccfloat_clear(x, ctx1);
        deccfloat_clear(y, ctx2);
        deccfloat_clear(z, ctx1);
        decfloat_clear(r, rctx);
        gr_ctx_clear(ctx1);
        gr_ctx_clear(ctx2);
        gr_ctx_clear(rctx);
    }

    TEST_FUNCTION_END(state);
}
