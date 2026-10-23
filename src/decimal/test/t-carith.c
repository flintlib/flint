/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "decimal.h"
#include "fmpq.h"
#include "acb.h"
#include "gr.h"

/* random finite complex number with moderate exponents */
static void
_crandtest_moderate(deccfloat_t x, flint_rand_t state, gr_ctx_t ctx)
{
    decfloat_ptr parts[2];
    int i;

    parts[0] = &x->re;
    parts[1] = &x->im;

    for (i = 0; i < 2; i++)
    {
        do
        {
            GR_MUST_SUCCEED(decfloat_randtest_special(parts[i], state, ctx));
        }
        while (!DECFLOAT_IS_FINITE(parts[i]));

        if (!DECFLOAT_IS_SPECIAL(parts[i]))
        {
            slong maxexp = 4 + 40 / DECIMAL_CTX_E(ctx);
            fmpz_set_si(&parts[i]->exp, (slong) n_randint(state, 2 * maxexp + 1) - maxexp);
        }
    }

    switch (n_randint(state, 6))
    {
        case 0: decfloat_zero(&x->im, ctx); break;
        case 1: decfloat_zero(&x->re, ctx); break;
        default: break;
    }
}

/* exact value as a pair of rationals */
static void
_get_fmpq2(fmpq_t re, fmpq_t im, const deccfloat_t x, gr_ctx_t ctx)
{
    GR_MUST_SUCCEED(decfloat_get_fmpq(re, &x->re, ctx));
    GR_MUST_SUCCEED(decfloat_get_fmpq(im, &x->im, ctx));
}

/* checks that y is the correct rounding of (re, im) */
static void
_check_fmpq2(const char * op, int status, const deccfloat_t y, const fmpq_t re, const fmpq_t im, const deccfloat_t x1, const deccfloat_t x2, gr_ctx_t ctx)
{
    deccfloat_t expected;
    int s2;

    deccfloat_init(expected, ctx);
    s2 = decfloat_set_round_fmpq_reference(&expected->re, re, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
    s2 |= decfloat_set_round_fmpq_reference(&expected->im, im, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND_IM(ctx), ctx);

    if (s2 != status || (status == GR_SUCCESS && deccfloat_equal(y, expected, ctx) != T_TRUE))
    {
        flint_printf("FAIL: %s (status %d, expected %d)\n", op, status, s2);
        gr_ctx_println(ctx);
        flint_printf("x1 = %{gr}\n", x1, ctx);
        flint_printf("x2 = %{gr}\n", x2, ctx);
        flint_printf("y = %{gr}\n", y, ctx);
        flint_printf("expected = %{gr}\n", expected, ctx);
        flint_abort();
    }

    deccfloat_clear(expected, ctx);
}

/* checks that the components of y are the correct roundings of the
   high-precision enclosure fa (which decides the rounding, or if it does
   not, y must at least lie within it) */
static void
_ccheck_acb(const char * op, int status, const deccfloat_t y, const acb_t fa, slong refprec, const deccfloat_t x, gr_ctx_t ctx)
{
    gr_ctx_t bctx;
    decball_t Y;
    decfloat_t lo, hi, r;
    slong prec = DECIMAL_CTX_PREC(ctx);
    int comp, decided = 1;

    if (!acb_is_finite(fa))
        return;

    _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), refprec, DECIMAL_RND_NEAR, 0);
    decball_init(Y, bctx);
    decfloat_init(lo, ctx);
    decfloat_init(hi, ctx);
    decfloat_init(r, bctx);

    /* first pass: does the enclosure decide both components? */
    for (comp = 0; comp < 2; comp++)
    {
        arb_srcptr part = comp ? acb_imagref(fa) : acb_realref(fa);
        int rnd = comp ? DECIMAL_CTX_RND_IM(ctx) : DECIMAL_CTX_RND(ctx);

        if (decball_set_arb(Y, part, bctx) != GR_SUCCESS)
        {
            flint_printf("FAIL: set_arb of reference\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(_decmag_get_decfloat(r, &Y->rad, bctx));
        GR_MUST_SUCCEED(decfloat_sub_round(lo, &Y->mid, r, refprec, DECIMAL_RND_FLOOR, bctx));
        GR_MUST_SUCCEED(decfloat_add_round(hi, &Y->mid, r, refprec, DECIMAL_RND_CEIL, bctx));
        GR_MUST_SUCCEED(decfloat_set_round(lo, lo, prec, rnd, ctx));
        GR_MUST_SUCCEED(decfloat_set_round(hi, hi, prec, rnd, ctx));

        if (decfloat_equal(lo, hi, ctx) != T_TRUE)
            decided = 0;
    }

    for (comp = 0; comp < 2; comp++)
    {
        arb_srcptr part = comp ? acb_imagref(fa) : acb_realref(fa);
        const decfloat_struct * yc = comp ? &y->im : &y->re;
        int rnd = comp ? DECIMAL_CTX_RND_IM(ctx) : DECIMAL_CTX_RND(ctx);

        GR_MUST_SUCCEED(decball_set_arb(Y, part, bctx));
        GR_MUST_SUCCEED(_decmag_get_decfloat(r, &Y->rad, bctx));
        GR_MUST_SUCCEED(decfloat_sub_round(lo, &Y->mid, r, refprec, DECIMAL_RND_FLOOR, bctx));
        GR_MUST_SUCCEED(decfloat_add_round(hi, &Y->mid, r, refprec, DECIMAL_RND_CEIL, bctx));
        GR_MUST_SUCCEED(decfloat_set_round(lo, lo, prec, rnd, ctx));
        GR_MUST_SUCCEED(decfloat_set_round(hi, hi, prec, rnd, ctx));

        if (decfloat_equal(lo, hi, ctx) == T_TRUE)
        {
            /* the rounding is decided: the function must agree (and must
               succeed if both components are decided) */
            if ((status == GR_SUCCESS && decfloat_equal(lo, yc, ctx) != T_TRUE) || (status != GR_SUCCESS && decided))
            {
                flint_printf("FAIL: %s (status %d)\n", op, status);
                gr_ctx_println(ctx);
                flint_printf("x = %{gr}\n", x, ctx);
                flint_printf("y = %{gr}\n", y, ctx);
                flint_printf("expected (component %d) = %s\n", comp, decfloat_get_str(lo, ctx));
                flint_printf("reference = %{acb}\n", fa);
                flint_abort();
            }
        }
        else if (status == GR_SUCCESS && (!_decfloat_is_finite(yc) ||
            _decfloat_cmp(yc, lo, ctx) < 0 || _decfloat_cmp(yc, hi, ctx) > 0))
        {
            flint_printf("FAIL: %s (undecided enclosure, status %d)\n", op, status);
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}\n", x, ctx);
            flint_printf("y = %{gr}\n", y, ctx);
            flint_printf("lo = %s, hi = %s\n", decfloat_get_str(lo, ctx), decfloat_get_str(hi, ctx));
            flint_printf("reference = %{acb}\n", fa);
            flint_abort();
        }
    }

    decball_clear(Y, bctx);
    decfloat_clear(lo, ctx);
    decfloat_clear(hi, ctx);
    decfloat_clear(r, bctx);
    gr_ctx_clear(bctx);
}

TEST_FUNCTION_START(deccfloat_arith, state)
{
    slong iter;

    /* exact operations against rationals */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        deccfloat_t x, y, z;
        fmpq_t a, b, c, d, re, im, t, u;
        int status, op, alias;

        gr_ctx_init_deccfloat_randtest(ctx, state, 30);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        if (DECIMAL_CTX_IS_EXACT(ctx) && n_randint(state, 2))
            decimal_ctx_set_prec(ctx, 1 + n_randint(state, 30));

        deccfloat_init(x, ctx);
        deccfloat_init(y, ctx);
        deccfloat_init(z, ctx);
        fmpq_init(a); fmpq_init(b); fmpq_init(c); fmpq_init(d);
        fmpq_init(re); fmpq_init(im); fmpq_init(t); fmpq_init(u);

        _crandtest_moderate(x, state, ctx);
        _crandtest_moderate(y, state, ctx);
        _get_fmpq2(a, b, x, ctx);
        _get_fmpq2(c, d, y, ctx);

        op = n_randint(state, 8);
        alias = n_randint(state, 3);

        switch (op)
        {
            case 0:
                fmpq_add(re, a, c);
                fmpq_add(im, b, d);
                if (alias == 1) { GR_MUST_SUCCEED(deccfloat_set_round(z, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, DECIMAL_RND_DOWN, ctx)); status = deccfloat_add(z, z, y, ctx); }
                else if (alias == 2) { GR_MUST_SUCCEED(deccfloat_set_round(z, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, DECIMAL_RND_DOWN, ctx)); status = deccfloat_add(z, x, z, ctx); }
                else status = deccfloat_add(z, x, y, ctx);
                _check_fmpq2("add", status, z, re, im, x, y, ctx);
                break;
            case 1:
                fmpq_sub(re, a, c);
                fmpq_sub(im, b, d);
                status = deccfloat_sub(z, x, y, ctx);
                _check_fmpq2("sub", status, z, re, im, x, y, ctx);
                break;
            case 2:
                /* (a + bi)(c + di) */
                fmpq_mul(re, a, c); fmpq_mul(t, b, d); fmpq_sub(re, re, t);
                fmpq_mul(im, a, d); fmpq_mul(t, b, c); fmpq_add(im, im, t);
                if (alias == 1) { GR_MUST_SUCCEED(deccfloat_set_round(z, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, DECIMAL_RND_DOWN, ctx)); status = deccfloat_mul(z, z, y, ctx); }
                else if (alias == 2) { GR_MUST_SUCCEED(deccfloat_set_round(z, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, DECIMAL_RND_DOWN, ctx)); status = deccfloat_mul(z, x, z, ctx); }
                else status = deccfloat_mul(z, x, y, ctx);
                _check_fmpq2("mul", status, z, re, im, x, y, ctx);
                break;
            case 3:
                fmpq_mul(re, a, a); fmpq_mul(t, b, b); fmpq_sub(re, re, t);
                fmpq_mul(im, a, b); fmpq_add(im, im, im);
                if (alias == 1) { GR_MUST_SUCCEED(deccfloat_set_round(z, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, DECIMAL_RND_DOWN, ctx)); status = deccfloat_sqr(z, z, ctx); }
                else status = deccfloat_sqr(z, x, ctx);
                _check_fmpq2("sqr", status, z, re, im, x, x, ctx);
                break;
            case 4:
            case 5:
                if (fmpq_is_zero(c) && fmpq_is_zero(d))
                    break;
                /* (a + bi) / (c + di) = ((ac + bd) + (bc - ad) i) / (c^2 + d^2) */
                fmpq_mul(u, c, c); fmpq_mul(t, d, d); fmpq_add(u, u, t);
                fmpq_mul(re, a, c); fmpq_mul(t, b, d); fmpq_add(re, re, t); fmpq_div(re, re, u);
                fmpq_mul(im, b, c); fmpq_mul(t, a, d); fmpq_sub(im, im, t); fmpq_div(im, im, u);
                if (alias == 1) { GR_MUST_SUCCEED(deccfloat_set_round(z, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, DECIMAL_RND_DOWN, ctx)); status = deccfloat_div(z, z, y, ctx); }
                else if (alias == 2) { GR_MUST_SUCCEED(deccfloat_set_round(z, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, DECIMAL_RND_DOWN, ctx)); status = deccfloat_div(z, x, z, ctx); }
                else status = deccfloat_div(z, x, y, ctx);
                _check_fmpq2("div", status, z, re, im, x, y, ctx);
                break;
            case 6:
                if (fmpq_is_zero(a) && fmpq_is_zero(b))
                    break;
                fmpq_mul(u, a, a); fmpq_mul(t, b, b); fmpq_add(u, u, t);
                fmpq_div(re, a, u);
                fmpq_div(im, b, u); fmpq_neg(im, im);
                status = deccfloat_inv(z, x, ctx);
                _check_fmpq2("inv", status, z, re, im, x, x, ctx);
                break;
            default:
                {
                    /* small integer powers */
                    slong n = (slong) n_randint(state, 13) - 6, i;
                    fmpq_t pr, pi;
                    fmpq_init(pr);
                    fmpq_init(pi);
                    fmpq_one(pr);
                    for (i = 0; i < FLINT_ABS(n); i++)
                    {
                        fmpq_mul(t, pr, a); fmpq_mul(u, pi, b); fmpq_sub(t, t, u);
                        fmpq_mul(u, pr, b); fmpq_mul(pi, pi, a); fmpq_add(pi, pi, u);
                        fmpq_set(pr, t);
                    }
                    if (n < 0)
                    {
                        if (fmpq_is_zero(pr) && fmpq_is_zero(pi))
                        {
                            fmpq_clear(pr);
                            fmpq_clear(pi);
                            break;
                        }
                        fmpq_mul(u, pr, pr); fmpq_mul(t, pi, pi); fmpq_add(u, u, t);
                        fmpq_div(re, pr, u);
                        fmpq_div(im, pi, u); fmpq_neg(im, im);
                    }
                    else
                    {
                        fmpq_set(re, pr);
                        fmpq_set(im, pi);
                    }
                    status = deccfloat_pow_si(z, x, n, ctx);
                    _check_fmpq2("pow_si", status, z, re, im, x, x, ctx);
                    fmpq_clear(pr);
                    fmpq_clear(pi);
                }
                break;
        }

        /* comparisons of absolute values */
        {
            int c1, c2, s1;
            fmpq_mul(t, a, a); fmpq_mul(u, b, b); fmpq_add(t, t, u);
            fmpq_mul(u, c, c); fmpq_mul(re, d, d); fmpq_add(u, u, re);
            c2 = fmpq_cmp(t, u);
            s1 = deccfloat_cmpabs(&c1, x, y, ctx);
            if (s1 != GR_SUCCESS || c1 != c2)
            {
                flint_printf("FAIL: cmpabs\n");
                gr_ctx_println(ctx);
                flint_printf("x = %{gr}\n", x, ctx);
                flint_printf("y = %{gr}\n", y, ctx);
                flint_printf("c1 = %d, c2 = %d\n", c1, c2);
                flint_abort();
            }
        }

        deccfloat_clear(x, ctx);
        deccfloat_clear(y, ctx);
        deccfloat_clear(z, ctx);
        fmpq_clear(a); fmpq_clear(b); fmpq_clear(c); fmpq_clear(d);
        fmpq_clear(re); fmpq_clear(im); fmpq_clear(t); fmpq_clear(u);
        gr_ctx_clear(ctx);
    }

    /* irrational operations against high-precision acb enclosures */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        deccfloat_t x, y, z;
        acb_t a, fa, b;
        slong prec, refprec, wpbits;
        int status, op;
        const char * name;

        gr_ctx_init_deccfloat_randtest(ctx, state, 20);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        if (DECIMAL_CTX_IS_EXACT(ctx))
            decimal_ctx_set_prec(ctx, 1 + n_randint(state, 20));
        prec = DECIMAL_CTX_PREC(ctx);

        deccfloat_init(x, ctx);
        deccfloat_init(y, ctx);
        deccfloat_init(z, ctx);
        acb_init(a);
        acb_init(fa);
        acb_init(b);

        _crandtest_moderate(x, state, ctx);
        _crandtest_moderate(y, state, ctx);

        /* exact squares now and then, to exercise the exact paths */
        if (n_randint(state, 4) == 0)
            GR_MUST_SUCCEED(deccfloat_sqr(x, x, ctx));

        refprec = 4 * prec + 100;
        wpbits = (slong) (3.33 * refprec) + 64;

        GR_MUST_SUCCEED(deccfloat_get_acb(a, x, wpbits, ctx));
        GR_MUST_SUCCEED(deccfloat_get_acb(b, y, wpbits, ctx));

        op = n_randint(state, 6);

        switch (op)
        {
            case 0:
                name = "sqrt";
                status = deccfloat_sqrt(z, x, ctx);
                acb_sqrt(fa, a, wpbits);
                break;
            case 1:
                name = "rsqrt";
                status = deccfloat_rsqrt(z, x, ctx);
                acb_rsqrt(fa, a, wpbits);
                break;
            case 2:
                name = "abs";
                status = deccfloat_abs(z, x, ctx);
                acb_abs(acb_realref(fa), a, wpbits);
                arb_zero(acb_imagref(fa));
                break;
            case 3:
                name = "sgn";
                status = deccfloat_sgn(z, x, ctx);
                acb_sgn(fa, a, wpbits);
                break;
            case 4:
                name = "arg";
                status = deccfloat_arg(z, x, ctx);
                acb_arg(acb_realref(fa), a, wpbits);
                arb_zero(acb_imagref(fa));
                break;
            default:
                name = "pow";
                /* small exponents so that the reference is accurate */
                if (n_randint(state, 2))
                    decfloat_zero(&y->im, ctx);
                if (n_randint(state, 2))
                    GR_MUST_SUCCEED(deccfloat_set_str(y, n_randint(state, 2) ? "0.5" : "-0.5", ctx));
                GR_MUST_SUCCEED(deccfloat_get_acb(b, y, wpbits, ctx));
                status = deccfloat_pow(z, x, y, ctx);
                acb_pow(fa, a, b, wpbits);
                break;
        }

        /* undefined values: any status */
        if (acb_is_finite(fa))
        {
            if (status == GR_SUCCESS)
                _ccheck_acb(name, status, z, fa, refprec, x, ctx);
            else if (status != GR_DOMAIN && status != GR_UNABLE)
            {
                flint_printf("FAIL: %s status %d\n", name, status);
                flint_abort();
            }

            /* a clean argument must succeed */
            if (status != GR_SUCCESS && op != 5 && deccfloat_is_zero(x, ctx) != T_TRUE && acb_rel_accuracy_bits(fa) > 3 * prec + 60)
            {
                flint_printf("FAIL: %s failed with status %d\n", name, status);
                gr_ctx_println(ctx);
                flint_printf("x = %{gr}\n", x, ctx);
                flint_printf("reference = %{acb}\n", fa);
                flint_abort();
            }
        }

        deccfloat_clear(x, ctx);
        deccfloat_clear(y, ctx);
        deccfloat_clear(z, ctx);
        acb_clear(a);
        acb_clear(fa);
        acb_clear(b);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
