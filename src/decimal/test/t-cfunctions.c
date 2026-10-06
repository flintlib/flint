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
#include "gr_special.h"

/* checks that the components of y are the correct roundings of the
   high-precision enclosure fa; returns 1 if the enclosure decided both
   components */
static int
_cfcheck_acb(slong which, slong kind, int status, const deccfloat_t y, const acb_t fa, slong refprec, const deccfloat_t x, gr_ctx_t ctx)
{
    gr_ctx_t bctx;
    decball_t Y;
    decfloat_t lo, hi, r;
    slong prec = DECIMAL_CTX_PREC(ctx);
    int comp, decided = 1;
    (void) kind;

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
                flint_printf("FAIL: function %wd (status %d)\n", which, status);
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
            flint_printf("FAIL: function %wd (undecided enclosure, status %d)\n", which, status);
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
    return decided;
}

/* random decfloat of moderate size: m 10^e with |m| < 10^6 and -8 <= e <= -2 */
static void
_randtest_moderate_part(decfloat_t x, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_t m, e;
    fmpz_init(m);
    fmpz_init(e);
    fmpz_set_ui(m, n_randint(state, 1000000));
    if (n_randint(state, 2)) fmpz_neg(m, m);
    fmpz_set_si(e, (slong) n_randint(state, 7) - 8);
    GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(x, m, e, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
    fmpz_clear(m);
    fmpz_clear(e);
}

TEST_FUNCTION_START(deccfloat_functions, state)
{
    slong iter;

    /* unary functions */
    {
        static const struct { int method; int near_one; int tiny; } tab[] = {
            {GR_METHOD_EXP, 0, 1}, {GR_METHOD_EXPM1, 0, 1}, {GR_METHOD_LOG, 1, 1},
            {GR_METHOD_LOG1P, 0, 1}, {GR_METHOD_LOG2, 1, 0}, {GR_METHOD_LOG10, 1, 0},
            {GR_METHOD_SIN, 0, 1}, {GR_METHOD_COS, 0, 1}, {GR_METHOD_TAN, 0, 1},
            {GR_METHOD_COT, 0, 1}, {GR_METHOD_SEC, 0, 1}, {GR_METHOD_CSC, 0, 1},
            {GR_METHOD_SINC, 0, 1}, {GR_METHOD_SIN_PI, 0, 0}, {GR_METHOD_COS_PI, 0, 1},
            {GR_METHOD_TAN_PI, 0, 0}, {GR_METHOD_COT_PI, 0, 0}, {GR_METHOD_SINC_PI, 0, 1},
            {GR_METHOD_SINH, 0, 1}, {GR_METHOD_COSH, 0, 1}, {GR_METHOD_TANH, 0, 1},
            {GR_METHOD_COTH, 0, 1}, {GR_METHOD_SECH, 0, 1}, {GR_METHOD_CSCH, 0, 1},
            {GR_METHOD_ASIN, 0, 1}, {GR_METHOD_ACOS, 0, 0}, {GR_METHOD_ATAN, 0, 1},
            {GR_METHOD_ACOT, 0, 0}, {GR_METHOD_ASINH, 0, 1}, {GR_METHOD_ACOSH, 0, 0},
            {GR_METHOD_ATANH, 0, 1}, {GR_METHOD_LAMBERTW, 0, 1}, {GR_METHOD_GAMMA, 0, 1},
            {GR_METHOD_RGAMMA, 0, 1}, {GR_METHOD_LGAMMA, 0, 0}, {GR_METHOD_DIGAMMA, 0, 1},
            {GR_METHOD_ZETA, 1, 1}, {GR_METHOD_ERF, 0, 0}, {GR_METHOD_ERFC, 0, 0},
            {GR_METHOD_ERFI, 0, 0}, {GR_METHOD_EXP_INTEGRAL_EI, 0, 0}, {GR_METHOD_SIN_INTEGRAL, 0, 1},
            {GR_METHOD_COS_INTEGRAL, 0, 0}, {GR_METHOD_SINH_INTEGRAL, 0, 1}, {GR_METHOD_DILOG, 0, 1},
            {GR_METHOD_AIRY_AI, 0, 0}, {GR_METHOD_AGM1, 0, 0}, {GR_METHOD_EXP_PI_I, 0, 0},
            {GR_METHOD_SQRT, 0, 0}, {GR_METHOD_RSQRT, 0, 0},
        };
        slong ntab = sizeof(tab) / sizeof(tab[0]);

        for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
        {
            gr_ctx_t ctx, bctx, actx;
            deccfloat_t x, y;
            deccball_t X, Y, Z;
            acb_t a, fa;
            fmpz_t m, t;
            slong prec, k, refprec, wpbits, which, kind;
            int status, ok;

            gr_ctx_init_deccfloat_randtest(ctx, state, 15);
            decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
            if (DECIMAL_CTX_IS_EXACT(ctx))
                decimal_ctx_set_prec(ctx, 1 + n_randint(state, 15));
            prec = DECIMAL_CTX_PREC(ctx);

            deccfloat_init(x, ctx);
            deccfloat_init(y, ctx);
            acb_init(a);
            acb_init(fa);
            fmpz_init(m);
            fmpz_init(t);

            which = n_randint(state, ntab);
            kind = n_randint(state, 5);
            k = 0;

            switch (kind)
            {
                case 0:
                    /* tiny: both parts +/- m 10^-k */
                    if (!tab[which].tiny)
                        kind = 4;
                    k = prec + 60 + n_randint(state, 100);
                    fmpz_set_ui(m, 1 + n_randint(state, 100000));
                    if (n_randint(state, 4) == 0) fmpz_one(m);
                    if (n_randint(state, 2)) fmpz_neg(m, m);
                    fmpz_set_si(t, -k);
                    GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(&x->re, m, t, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                    if (n_randint(state, 3) != 0)
                    {
                        fmpz_set_ui(m, 1 + n_randint(state, 100000));
                        if (n_randint(state, 4) == 0) fmpz_one(m);
                        fmpz_set_si(t, -k - n_randint(state, 3));
                    }
                    if (n_randint(state, 2)) fmpz_neg(m, m);
                    GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(&x->im, m, t, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                    break;
                case 1:
                    /* real rational points */
                    {
                        fmpq_t q;
                        fmpq_init(q);
                        if (n_randint(state, 2))
                            fmpq_set_si(q, (slong) n_randint(state, 49) - 24, 4);
                        else
                            fmpq_set_si(q, (slong) n_randint(state, 21) - 10, 1 + n_randint(state, 2));
                        GR_MUST_SUCCEED(decfloat_set_round_fmpq(&x->re, q, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                        decfloat_zero(&x->im, ctx);
                        fmpq_clear(q);
                    }
                    break;
                case 2:
                    /* imaginary */
                    _randtest_moderate_part(&x->im, state, ctx);
                    decfloat_zero(&x->re, ctx);
                    if (n_randint(state, 4) == 0)
                    {
                        /* tiny imaginary */
                        k = prec + 60 + n_randint(state, 100);
                        fmpz_set_ui(m, 1 + n_randint(state, 100000));
                        fmpz_set_si(t, -k);
                        GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(&x->im, m, t, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                    }
                    break;
                case 3:
                    /* real random */
                    _randtest_moderate_part(&x->re, state, ctx);
                    decfloat_zero(&x->im, ctx);
                    break;
                default:
                    /* moderate random */
                    _randtest_moderate_part(&x->re, state, ctx);
                    _randtest_moderate_part(&x->im, state, ctx);
                    break;
            }

            if (tab[which].near_one)
            {
                decfloat_t one;
                decfloat_init(one, ctx);
                GR_MUST_SUCCEED(decfloat_one(one, ctx));
                GR_MUST_SUCCEED(decfloat_add_round(&x->re, &x->re, one, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                decfloat_clear(one, ctx);
            }

            if (!_deccfloat_is_finite(x))
                goto cleanup;

            status = ((gr_method_unary_op) ctx->methods[tab[which].method])(y, x, ctx);

            refprec = 4 * k + 4 * prec + 100;
            wpbits = (slong) (3.33 * refprec) + 64;
            gr_ctx_init_complex_acb(actx, wpbits);

            GR_MUST_SUCCEED(deccfloat_get_acb(a, x, wpbits, ctx));

            ok = (((gr_method_unary_op) actx->methods[tab[which].method])(fa, a, actx) == GR_SUCCESS) && acb_is_finite(fa);

            if (ok)
            {
                int decided = _cfcheck_acb(which, kind, status, y, fa, refprec, x, ctx);

                /* tiny arguments must always be handled */
                if (kind == 0 && status != GR_SUCCESS && decided)
                {
                    flint_printf("FAIL: tiny argument not handled, function %wd, status = %d\n", which, status);
                    gr_ctx_println(ctx);
                    flint_printf("x = %{gr}\n", x, ctx);
                    flint_printf("reference = %{acb}\n", fa);
                    flint_abort();
                }

                /* the ball function must contain the high-precision value */
                _gr_ctx_init_decimal(bctx, DECIMAL_CTX_CBALL, DECIMAL_CTX_E(ctx), prec + 40, DECIMAL_RND_NEAR, 0);
                deccball_init(X, bctx);
                deccball_init(Y, bctx);
                deccball_init(Z, bctx);
                GR_MUST_SUCCEED(deccball_set_deccfloat(X, x, bctx));
                if (((gr_method_unary_op) bctx->methods[tab[which].method])(Y, X, bctx) == GR_SUCCESS)
                {
                    GR_MUST_SUCCEED(deccball_set_acb(Z, fa, bctx));
                    if (!_deccball_overlaps(Y, Z, bctx))
                    {
                        flint_printf("FAIL: ball function %wd, kind %wd\n", which, kind);
                        gr_ctx_println(bctx);
                        flint_printf("X = %{gr}\n", X, bctx);
                        flint_printf("Y = %{gr}\n", Y, bctx);
                        flint_printf("Z = %{gr}\n", Z, bctx);
                        flint_abort();
                    }
                }
                deccball_clear(X, bctx);
                deccball_clear(Y, bctx);
                deccball_clear(Z, bctx);
                gr_ctx_clear(bctx);
            }

            gr_ctx_clear(actx);

cleanup:
            deccfloat_clear(x, ctx);
            deccfloat_clear(y, ctx);
            acb_clear(a);
            acb_clear(fa);
            fmpz_clear(m);
            fmpz_clear(t);
            gr_ctx_clear(ctx);
        }
    }

    /* a few binary functions */
    {
        static const struct { int method; } tab[] = {
            {GR_METHOD_AGM}, {GR_METHOD_BESSEL_J}, {GR_METHOD_BESSEL_K},
            {GR_METHOD_HURWITZ_ZETA}, {GR_METHOD_RISING}, {GR_METHOD_POW},
        };
        slong ntab = sizeof(tab) / sizeof(tab[0]);

        for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
        {
            gr_ctx_t ctx, bctx, actx;
            deccfloat_t x, y, z;
            deccball_t X, Y, Z, W;
            acb_t a, b, fa;
            slong prec, refprec, wpbits, which;
            int status, ok;

            gr_ctx_init_deccfloat_randtest(ctx, state, 15);
            decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
            if (DECIMAL_CTX_IS_EXACT(ctx))
                decimal_ctx_set_prec(ctx, 1 + n_randint(state, 15));
            prec = DECIMAL_CTX_PREC(ctx);

            deccfloat_init(x, ctx);
            deccfloat_init(y, ctx);
            deccfloat_init(z, ctx);
            acb_init(a);
            acb_init(b);
            acb_init(fa);

            which = n_randint(state, ntab);

            _randtest_moderate_part(&x->re, state, ctx);
            _randtest_moderate_part(&x->im, state, ctx);
            _randtest_moderate_part(&y->re, state, ctx);
            _randtest_moderate_part(&y->im, state, ctx);
            if (n_randint(state, 4) == 0) decfloat_zero(&x->im, ctx);
            if (n_randint(state, 4) == 0) decfloat_zero(&y->im, ctx);
            if (n_randint(state, 8) == 0) decfloat_zero(&x->re, ctx);
            if (n_randint(state, 8) == 0) decfloat_zero(&y->re, ctx);

            /* avoid exact cases (integer orders, zero arguments), where
               Ziv's loop runs to its limit without deciding */
            if (deccfloat_is_zero(x, ctx) == T_TRUE || deccfloat_is_zero(y, ctx) == T_TRUE)
                goto cleanup2;
            if (which == 3)
            {
                /* hurwitz zeta: moderate s away from the integers, where
                   the acb implementation is fast */
                if (_decfloat_is_int(&x->re, ctx))
                    goto cleanup2;
                GR_MUST_SUCCEED(decfloat_mul_10exp_si_round(&x->re, &x->re, -3, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                GR_MUST_SUCCEED(decfloat_mul_10exp_si_round(&x->im, &x->im, -3, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
            }
            if (which == 4 && _decfloat_is_int(&y->re, ctx) && DECFLOAT_IS_ZERO(&y->im))
                goto cleanup2;

            if (!_deccfloat_is_finite(x) || !_deccfloat_is_finite(y))
                goto cleanup2;

            status = ((gr_method_binary_op) ctx->methods[tab[which].method])(z, x, y, ctx);

            refprec = 4 * prec + 100;
            wpbits = (slong) (3.33 * refprec) + 64;
            gr_ctx_init_complex_acb(actx, wpbits);

            GR_MUST_SUCCEED(deccfloat_get_acb(a, x, wpbits, ctx));
            GR_MUST_SUCCEED(deccfloat_get_acb(b, y, wpbits, ctx));

            ok = (((gr_method_binary_op) actx->methods[tab[which].method])(fa, a, b, actx) == GR_SUCCESS) && acb_is_finite(fa);

            if (ok)
            {
                _cfcheck_acb(which, 10, status, z, fa, refprec, x, ctx);

                _gr_ctx_init_decimal(bctx, DECIMAL_CTX_CBALL, DECIMAL_CTX_E(ctx), prec + 40, DECIMAL_RND_NEAR, 0);
                deccball_init(X, bctx);
                deccball_init(Y, bctx);
                deccball_init(Z, bctx);
                deccball_init(W, bctx);
                GR_MUST_SUCCEED(deccball_set_deccfloat(X, x, bctx));
                GR_MUST_SUCCEED(deccball_set_deccfloat(Y, y, bctx));
                if (((gr_method_binary_op) bctx->methods[tab[which].method])(Z, X, Y, bctx) == GR_SUCCESS)
                {
                    GR_MUST_SUCCEED(deccball_set_acb(W, fa, bctx));
                    if (!_deccball_overlaps(Z, W, bctx))
                    {
                        flint_printf("FAIL: binary ball function %wd\n", which);
                        gr_ctx_println(bctx);
                        flint_printf("X = %{gr}\n", X, bctx);
                        flint_printf("Y = %{gr}\n", Y, bctx);
                        flint_printf("Z = %{gr}\n", Z, bctx);
                        flint_printf("W = %{gr}\n", W, bctx);
                        flint_abort();
                    }
                }
                deccball_clear(X, bctx);
                deccball_clear(Y, bctx);
                deccball_clear(Z, bctx);
                deccball_clear(W, bctx);
                gr_ctx_clear(bctx);
            }

            gr_ctx_clear(actx);

cleanup2:
            deccfloat_clear(x, ctx);
            deccfloat_clear(y, ctx);
            deccfloat_clear(z, ctx);
            acb_clear(a);
            acb_clear(b);
            acb_clear(fa);
            gr_ctx_clear(ctx);
        }
    }

    /* regression: the tiny-argument tail bound must include the leading
       tail term (sinc(z) = 1 - z^2/6 + ..., with a small limb radix) */
    {
        gr_ctx_t ctx;
        gr_ptr x, y, z;

        _gr_ctx_init_decimal(ctx, DECIMAL_CTX_CFLOAT, 1, 11, DECIMAL_RND_CEIL, 0);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);
        GR_MUST_SUCCEED(gr_set_str(x, "-0.00036857*I", ctx));
        GR_MUST_SUCCEED(gr_set_str(z, "1.0000000227", ctx));
        if (gr_sinc(y, x, ctx) != GR_SUCCESS || gr_equal(y, z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: sinc (tiny argument)\n");
            flint_printf("y = %{gr}\n", y, ctx);
            flint_abort();
        }
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
