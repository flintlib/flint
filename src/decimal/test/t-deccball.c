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
#include "fmpq.h"
#include "acb.h"
#include "gr.h"
#include "gr_poly.h"

/* random ball with moderate exponents; stores a random rational point of
   the ball in q */
static void
_crandom_ball_point(decball_t x, fmpq_t q, flint_rand_t state, gr_ctx_t ctx)
{
    fmpq_t m, r, t;

    fmpq_init(m);
    fmpq_init(r);
    fmpq_init(t);

    do
    {
        GR_MUST_SUCCEED(decfloat_randtest_special(&x->mid, state, ctx));
    }
    while (!DECFLOAT_IS_FINITE(&x->mid));

    if (!DECFLOAT_IS_SPECIAL(&x->mid))
        fmpz_set_si(&x->mid.exp, (slong) n_randint(state, 5) - 2);

    switch (n_randint(state, 4))
    {
        case 0:
            _decmag_zero(&x->rad, ctx);
            break;
        case 1:
            {
                slong E;
                if (!DECFLOAT_IS_SPECIAL(&x->mid) && decfloat_get_sci_exp_si(&E, &x->mid, ctx))
                    _decmag_set_ui_10exp_si(&x->rad, 1 + n_randint(state, DECIMAL_CTX_RAD_POW(ctx) - 1),
                        E - DECIMAL_CTX_PREC(ctx) - DECIMAL_CTX_RAD_PREC(ctx) + n_randint(state, 4), ctx);
                else
                    _decmag_zero(&x->rad, ctx);
            }
            break;
        default:
            _decmag_set_ui_10exp_si(&x->rad, 1 + n_randint(state, DECIMAL_CTX_RAD_POW(ctx) - 1),
                (slong) n_randint(state, 21) - 15, ctx);
            break;
    }

    GR_MUST_SUCCEED(decfloat_get_fmpq(m, &x->mid, ctx));
    _decmag_get_fmpq(r, &x->rad, ctx);

    switch (n_randint(state, 6))
    {
        case 0: fmpq_set_si(t, 1, 1); break;
        case 1: fmpq_set_si(t, -1, 1); break;
        case 2: fmpq_zero(t); break;
        default:
            {
                fmpz_t a, b;
                fmpz_init(a);
                fmpz_init(b);
                fmpz_randtest(a, state, 30);
                fmpz_randtest_not_zero(b, state, 30);
                fmpz_abs(b, b);
                fmpz_fdiv_r(a, a, b);
                fmpq_set_fmpz_frac(t, a, b);
                if (n_randint(state, 2)) fmpq_neg(t, t);
                fmpz_clear(a);
                fmpz_clear(b);
            }
            break;
    }

    fmpq_mul(t, t, r);
    fmpq_add(q, m, t);

    fmpq_clear(m);
    fmpq_clear(r);
    fmpq_clear(t);
}

static void
_random_cball_point(deccball_t x, fmpq_t re, fmpq_t im, flint_rand_t state, gr_ctx_t ctx)
{
    _crandom_ball_point(&x->re, re, state, ctx);
    _crandom_ball_point(&x->im, im, state, ctx);

    switch (n_randint(state, 6))
    {
        case 0: decball_zero(&x->im, ctx); fmpq_zero(im); break;
        case 1: decball_zero(&x->re, ctx); fmpq_zero(re); break;
        default: break;
    }
}

static int
_contains_fmpq2(const deccball_t x, const fmpq_t re, const fmpq_t im, gr_ctx_t ctx)
{
    return _decball_contains_fmpq(&x->re, re, ctx) && _decball_contains_fmpq(&x->im, im, ctx);
}

/* parse s in ctx and check that the result prints as expected */
static void
_cbcheck_parse(const char * s, const char * expected, gr_ctx_t ctx)
{
    deccball_t x;
    char * t;
    int status;

    deccball_init(x, ctx);
    status = gr_set_str(x, s, ctx);
    t = (status == GR_SUCCESS) ? deccball_get_str(x, ctx) : NULL;

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
    deccball_clear(x, ctx);
}

TEST_FUNCTION_START(deccball, state)
{
    slong iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_ctx_init_deccball_randtest(ctx, state, 40);
        gr_test_ring(ctx, 30, 0);
        gr_ctx_clear(ctx);
    }

    {
        gr_ctx_t ctx;
        gr_ctx_init_deccball(ctx, 10, 0);
        gr_test_ring(ctx, 100 * flint_test_multiplier(), 0);
        gr_ctx_clear(ctx);
        gr_ctx_init_deccball(ctx, 30, DECIMAL_SLOPPY_RADIUS);
        gr_test_ring(ctx, 100 * flint_test_multiplier(), 0);
        gr_ctx_clear(ctx);
    }

    /* strings */
    {
        gr_ctx_t ctx;

        gr_ctx_init_deccball(ctx, 20, 0);
        decimal_ctx_set_rad_prec(ctx, 3);

        _cbcheck_parse("1.5", "1.5", ctx);
        _cbcheck_parse("[1.5 +/- 0.01]", "[1.5 +/- 0.01]", ctx);
        _cbcheck_parse("2*I", "2*I", ctx);
        _cbcheck_parse("(1 + 2*I)", "(1 + 2*I)", ctx);
        _cbcheck_parse("([1 +/- 0.1] + [2 +/- 0.2]*I)", "([1 +/- 0.1] + [2 +/- 0.2]*I)", ctx);
        _cbcheck_parse("([1 +/- 0.1] - 2*I)", "([1 +/- 0.1] - 2*I)", ctx);
        _cbcheck_parse("[1 +/- 0.1]*I", "[1 +/- 0.1]*I", ctx);
        _cbcheck_parse("[1 + 2*I +/- 0.1]", "([1 +/- 0.1] + 2*I)", ctx);
        _cbcheck_parse("[1 + 2*I +/- (0.1 + 0.2*I)]", "([1 +/- 0.1] + [2 +/- 0.2]*I)", ctx);
        _cbcheck_parse("(1 + 2*I) * (3 - 4*I)", "(11 + 2*I)", ctx);
        _cbcheck_parse("(1 + 2*I) / (3 - 4*I)", "(-0.2 + 0.4*I)", ctx);
        _cbcheck_parse("1/3 + I/3", "([0.33333333333333333333 +/- 3.34e-21] + [0.33333333333333333333 +/- 3.34e-21]*I)", ctx);
        _cbcheck_parse("sqrt(-4)", "2*I", ctx);
        _cbcheck_parse("(3+4*I)^(1/2)", "(2 + 1*I)", ctx);

        gr_ctx_clear(ctx);
    }

    /* string roundtrips, also in polynomials */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, pctx;
        deccball_t x, y;
        gr_poly_t f, g;
        char * s;
        int status;
        slong i;

        gr_ctx_init_deccball_randtest(ctx, state, 30);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        gr_ctx_init_gr_poly(pctx, ctx);

        deccball_init(x, ctx);
        deccball_init(y, ctx);

        GR_MUST_SUCCEED(deccball_randtest(x, state, ctx));
        s = deccball_get_str(x, ctx);
        status = gr_set_str(y, s, ctx);

        if (status != GR_SUCCESS || !_deccball_contains(y, x, ctx))
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
        status |= gr_get_str(&s, f, pctx);
        status |= gr_set_str(g, s, pctx);

        if (status != GR_SUCCESS || f->length != g->length)
        {
            flint_printf("FAIL: polynomial roundtrip\n");
            gr_ctx_println(ctx);
            flint_printf("f = %s\n", s);
            flint_printf("g = %{gr}  (%d)\n", g, pctx, status);
            flint_abort();
        }

        for (i = 0; i < f->length; i++)
        {
            if (!_deccball_contains(gr_poly_coeff_ptr(g, i, ctx), gr_poly_coeff_ptr(f, i, ctx), ctx))
            {
                flint_printf("FAIL: polynomial roundtrip (containment)\n");
                gr_ctx_println(ctx);
                flint_printf("f = %s\n", s);
                flint_printf("g = %{gr}\n", g, pctx);
                flint_abort();
            }
        }

        flint_free(s);
        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        deccball_clear(x, ctx);
        deccball_clear(y, ctx);
        gr_ctx_clear(pctx);
        gr_ctx_clear(ctx);
    }

    /* containment of exact results */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        deccball_t X, Y, Z;
        fmpq_t a, b, c, d, re, im, t, u;
        int status, op;
        const char * name;

        gr_ctx_init_deccball_randtest(ctx, state, 30);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);

        deccball_init(X, ctx);
        deccball_init(Y, ctx);
        deccball_init(Z, ctx);
        fmpq_init(a); fmpq_init(b); fmpq_init(c); fmpq_init(d);
        fmpq_init(re); fmpq_init(im); fmpq_init(t); fmpq_init(u);

        _random_cball_point(X, a, b, state, ctx);
        _random_cball_point(Y, c, d, state, ctx);

        op = n_randint(state, 7);

        switch (op)
        {
            case 0:
                name = "add";
                fmpq_add(re, a, c); fmpq_add(im, b, d);
                status = deccball_add(Z, X, Y, ctx);
                break;
            case 1:
                name = "sub";
                fmpq_sub(re, a, c); fmpq_sub(im, b, d);
                status = deccball_sub(Z, X, Y, ctx);
                break;
            case 2:
                name = "mul";
                fmpq_mul(re, a, c); fmpq_mul(t, b, d); fmpq_sub(re, re, t);
                fmpq_mul(im, a, d); fmpq_mul(t, b, c); fmpq_add(im, im, t);
                status = deccball_mul(Z, X, Y, ctx);
                break;
            case 3:
                name = "sqr";
                fmpq_mul(re, a, a); fmpq_mul(t, b, b); fmpq_sub(re, re, t);
                fmpq_mul(im, a, b); fmpq_add(im, im, im);
                status = deccball_sqr(Z, X, ctx);
                break;
            case 4:
                name = "div";
                fmpq_mul(u, c, c); fmpq_mul(t, d, d); fmpq_add(u, u, t);
                if (fmpq_is_zero(u))
                    goto next;
                fmpq_mul(re, a, c); fmpq_mul(t, b, d); fmpq_add(re, re, t); fmpq_div(re, re, u);
                fmpq_mul(im, b, c); fmpq_mul(t, a, d); fmpq_sub(im, im, t); fmpq_div(im, im, u);
                status = deccball_div(Z, X, Y, ctx);
                break;
            case 5:
                name = "inv";
                fmpq_mul(u, a, a); fmpq_mul(t, b, b); fmpq_add(u, u, t);
                if (fmpq_is_zero(u))
                    goto next;
                fmpq_div(re, a, u);
                fmpq_div(im, b, u); fmpq_neg(im, im);
                status = deccball_inv(Z, X, ctx);
                break;
            default:
                name = "pow_ui";
                {
                    ulong n = n_randint(state, 6), i;
                    fmpq_one(re);
                    fmpq_zero(im);
                    for (i = 0; i < n; i++)
                    {
                        fmpq_mul(t, re, a); fmpq_mul(u, im, b); fmpq_sub(t, t, u);
                        fmpq_mul(u, re, b); fmpq_mul(im, im, a); fmpq_add(im, im, u);
                        fmpq_set(re, t);
                    }
                    status = deccball_pow_ui(Z, X, n, ctx);
                }
                break;
        }

        if (status == GR_SUCCESS && !_contains_fmpq2(Z, re, im, ctx))
        {
            flint_printf("FAIL: containment (%s)\n", name);
            gr_ctx_println(ctx);
            flint_printf("X = %{gr}\n", X, ctx);
            flint_printf("Y = %{gr}\n", Y, ctx);
            flint_printf("Z = %{gr}\n", Z, ctx);
            flint_printf("re = %{fmpq}, im = %{fmpq}\n", re, im);
            flint_abort();
        }

        if (status != GR_SUCCESS && (op <= 3 || op == 6))
        {
            flint_printf("FAIL: status %d (%s)\n", status, name);
            gr_ctx_println(ctx);
            flint_printf("X = %{gr}\n", X, ctx);
            flint_printf("Y = %{gr}\n", Y, ctx);
            flint_abort();
        }

next:
        deccball_clear(X, ctx);
        deccball_clear(Y, ctx);
        deccball_clear(Z, ctx);
        fmpq_clear(a); fmpq_clear(b); fmpq_clear(c); fmpq_clear(d);
        fmpq_clear(re); fmpq_clear(im); fmpq_clear(t); fmpq_clear(u);
        gr_ctx_clear(ctx);
    }

    /* containment of irrational results (compared with acb) */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        deccball_t X, Y, Z, W;
        acb_t a, b, fa;
        fmpq_t qa, qb, qc, qd;
        slong wp;
        int status, op;
        const char * name;

        gr_ctx_init_deccball_randtest(ctx, state, 30);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);

        deccball_init(X, ctx);
        deccball_init(Y, ctx);
        deccball_init(Z, ctx);
        deccball_init(W, ctx);
        acb_init(a);
        acb_init(b);
        acb_init(fa);
        fmpq_init(qa); fmpq_init(qb); fmpq_init(qc); fmpq_init(qd);

        _random_cball_point(X, qa, qb, state, ctx);
        _random_cball_point(Y, qc, qd, state, ctx);

        /* the point */
        wp = _decimal_digits_to_bits(DECIMAL_CTX_PREC(ctx)) + 200;
        arb_set_fmpq(acb_realref(a), qa, wp);
        arb_set_fmpq(acb_imagref(a), qb, wp);
        arb_set_fmpq(acb_realref(b), qc, wp);
        arb_set_fmpq(acb_imagref(b), qd, wp);

        op = n_randint(state, 7);

        switch (op)
        {
            case 0: name = "sqrt"; status = deccball_sqrt(Z, X, ctx); acb_sqrt(fa, a, wp); break;
            case 1: name = "rsqrt"; status = deccball_rsqrt(Z, X, ctx); acb_rsqrt(fa, a, wp); break;
            case 2: name = "abs"; status = deccball_abs(Z, X, ctx); acb_abs(acb_realref(fa), a, wp); arb_zero(acb_imagref(fa)); break;
            case 3: name = "sgn"; status = deccball_sgn(Z, X, ctx); acb_sgn(fa, a, wp); break;
            case 4: name = "arg"; status = deccball_arg(Z, X, ctx); acb_arg(acb_realref(fa), a, wp); arb_zero(acb_imagref(fa)); break;
            case 5: name = "pow"; status = deccball_pow(Z, X, Y, ctx); acb_pow(fa, a, b, wp); break;
            default: name = "exp"; status = deccball_exp(Z, X, ctx); acb_exp(fa, a, wp); break;
        }

        if (status == GR_SUCCESS && acb_is_finite(fa))
        {
            GR_MUST_SUCCEED(deccball_set_acb(W, fa, ctx));

            if (!_deccball_overlaps(Z, W, ctx))
            {
                flint_printf("FAIL: containment (%s)\n", name);
                gr_ctx_println(ctx);
                flint_printf("X = %{gr}\n", X, ctx);
                flint_printf("Y = %{gr}\n", Y, ctx);
                flint_printf("Z = %{gr}\n", Z, ctx);
                flint_printf("W = %{gr}\n", W, ctx);
                flint_abort();
            }
        }

        deccball_clear(X, ctx);
        deccball_clear(Y, ctx);
        deccball_clear(Z, ctx);
        deccball_clear(W, ctx);
        acb_clear(a);
        acb_clear(b);
        acb_clear(fa);
        fmpq_clear(qa); fmpq_clear(qb); fmpq_clear(qc); fmpq_clear(qd);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
