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
#include "gr_special.h"
#include "fmpq.h"
#include "arb.h"
#include "gr.h"

TEST_FUNCTION_START(decimal_via_arb, state)
{
    slong iter;

    /* ball <-> arb conversions are rigorous */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decball_t X, Y;
        arb_t a, b;
        fmpq_t q;
        slong pb;

        gr_ctx_init_decball_randtest(ctx, state, 40);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);

        decball_init(X, ctx);
        decball_init(Y, ctx);
        arb_init(a);
        arb_init(b);
        fmpq_init(q);

        GR_MUST_SUCCEED(decball_randtest(X, state, ctx));
        if (_decball_is_finite(X, ctx) && fmpz_bits(&X->mid.exp) < 8 && fmpz_bits(&X->rad.exp) < 8)
        {
            pb = 2 + n_randint(state, 300);
            GR_MUST_SUCCEED(decball_get_arb(a, X, pb, ctx));

            /* the midpoint (a rational) must be contained in a */
            GR_MUST_SUCCEED(decfloat_get_fmpq(q, &X->mid, ctx));
            if (!arb_contains_fmpq(a, q))
            {
                flint_printf("FAIL: get_arb\n");
                flint_printf("X = %{gr}, a = %{arb}\n", X, ctx, a);
                flint_abort();
            }

            /* converting back must contain the original midpoint */
            GR_MUST_SUCCEED(decball_set_arb(Y, a, ctx));
            if (!_decball_contains_fmpq(Y, q, ctx))
            {
                flint_printf("FAIL: set_arb\n");
                flint_printf("X = %{gr}, a = %{arb}, Y = %{gr}\n", X, ctx, a, Y, ctx);
                flint_abort();
            }
        }

        /* arb -> decball -> contains random point of the arb ball */
        arb_randtest(a, state, 1 + n_randint(state, 200), 6);
        if (arb_is_finite(a))
        {
            arf_t p;
            arf_init(p);
            arb_get_mid_arb(b, a);
            arf_set_mag(p, arb_radref(a));
            if (n_randint(state, 2)) arf_neg(p, p);
            arf_add(p, arb_midref(a), p, ARF_PREC_EXACT, ARF_RND_DOWN);
            /* p is an endpoint */
            GR_MUST_SUCCEED(decball_set_arb(Y, a, ctx));
            {
                decfloat_t f;
                decfloat_init(f, ctx);
                if (decfloat_set_round_arf(f, p, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx) == GR_SUCCESS)
                {
                    if (!_decball_contains_decfloat(Y, f, ctx))
                    {
                        flint_printf("FAIL: set_arb endpoint\n");
                        flint_printf("a = %{arb}, Y = %{gr}, f = %{gr}\n", a, Y, ctx, f, ctx);
                        flint_abort();
                    }
                }
                decfloat_clear(f, ctx);
            }
            arf_clear(p);
        }

        decball_clear(X, ctx);
        decball_clear(Y, ctx);
        arb_clear(a);
        arb_clear(b);
        fmpq_clear(q);
        gr_ctx_clear(ctx);
    }

    /* transcendental functions on balls: containment via high-precision arb */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decball_t X, Y;
        arb_t a, b;
        int which, status;

        gr_ctx_init_decball_randtest(ctx, state, 30);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);

        decball_init(X, ctx);
        decball_init(Y, ctx);
        arb_init(a);
        arb_init(b);

        GR_MUST_SUCCEED(decball_randtest(X, state, ctx));
        if (!DECFLOAT_IS_SPECIAL(&X->mid))
            fmpz_set_si(&X->mid.exp, (slong) n_randint(state, 5) - 3);
        if (!DECMAG_IS_SPECIAL(&X->rad) && fmpz_bits(&X->rad.exp) > 6)
            _decmag_zero(&X->rad, ctx);

        which = n_randint(state, 8);

        switch (which)
        {
            case 0: status = decball_exp(Y, X, ctx); break;
            case 1: status = decball_log(Y, X, ctx); break;
            case 2: status = decball_sin(Y, X, ctx); break;
            case 3: status = decball_cos(Y, X, ctx); break;
            case 4: status = decball_atan(Y, X, ctx); break;
            case 5: status = decball_sqrt(Y, X, ctx); break;
            case 6: status = decball_expm1(Y, X, ctx); break;
            default: status = decball_tanh(Y, X, ctx); break;
        }

        if (status == GR_SUCCESS && _decball_is_finite(X, ctx) && _decball_is_finite(Y, ctx))
        {
            /* a = f(X) at high precision must overlap with Y (since Y contains f(x) for all x in X,
               and the arb ball also does) -- in fact Y must contain f(mid) */
            decball_t M;
            decball_init(M, ctx);
            GR_MUST_SUCCEED(decfloat_set_round(&M->mid, &X->mid, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
            _decmag_zero(&M->rad, ctx);
            GR_MUST_SUCCEED(decball_get_arb(a, M, 400, ctx));

            switch (which)
            {
                case 0: arb_exp(b, a, 400); break;
                case 1: arb_log(b, a, 400); break;
                case 2: arb_sin(b, a, 400); break;
                case 3: arb_cos(b, a, 400); break;
                case 4: arb_atan(b, a, 400); break;
                case 5: arb_sqrt(b, a, 400); break;
                case 6: arb_expm1(b, a, 400); break;
                default: arb_tanh(b, a, 400); break;
            }

            GR_MUST_SUCCEED(decball_get_arb(a, Y, 400, ctx));

            if (arb_is_finite(b) && !arb_overlaps(a, b))
            {
                flint_printf("FAIL: transcendental containment, which = %d\n", which);
                gr_ctx_println(ctx);
                flint_printf("X = %{gr}\n", X, ctx);
                flint_printf("Y = %{gr}\n", Y, ctx);
                flint_printf("b = %{arb}\n", b);
                flint_abort();
            }

            decball_clear(M, ctx);
        }

        decball_clear(X, ctx);
        decball_clear(Y, ctx);
        arb_clear(a);
        arb_clear(b);
        gr_ctx_clear(ctx);
    }

    /* correctly rounded float functions: check against the reference rounding
       of a tight rational enclosure */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, bctx;
        decfloat_t x, y, lo, hi;
        decball_t B;
        arb_t a;
        int which, status;
        slong prec;

        gr_ctx_init_decfloat_randtest(ctx, state, 30);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        if (DECIMAL_CTX_IS_EXACT(ctx))
            decimal_ctx_set_prec(ctx, 1 + n_randint(state, 30));
        prec = DECIMAL_CTX_PREC(ctx);

        _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec + 40, DECIMAL_RND_NEAR, 0);

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(lo, ctx);
        decfloat_init(hi, ctx);
        decball_init(B, bctx);
        arb_init(a);

        GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx));
        if (!DECFLOAT_IS_SPECIAL(x))
            fmpz_set_si(&x->exp, (slong) n_randint(state, 5) - 3);

        which = n_randint(state, 5);

        switch (which)
        {
            case 0: status = decfloat_exp(y, x, ctx); break;
            case 1: status = decfloat_log(y, x, ctx); break;
            case 2: status = decfloat_sin(y, x, ctx); break;
            case 3: status = decfloat_cos(y, x, ctx); break;
            default: status = decfloat_atan(y, x, ctx); break;
        }

        if (status == GR_SUCCESS && DECFLOAT_IS_FINITE(x) && DECFLOAT_IS_FINITE(y))
        {
            GR_MUST_SUCCEED(decfloat_set_round(&B->mid, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, bctx));
            _decmag_zero(&B->rad, bctx);
            GR_MUST_SUCCEED(decball_get_arb(a, B, 4 * prec + 200, bctx));

            switch (which)
            {
                case 0: arb_exp(a, a, 4 * prec + 200); break;
                case 1: arb_log(a, a, 4 * prec + 200); break;
                case 2: arb_sin(a, a, 4 * prec + 200); break;
                case 3: arb_cos(a, a, 4 * prec + 200); break;
                default: arb_atan(a, a, 4 * prec + 200); break;
            }

            if (arb_is_finite(a))
            {
                /* y must equal the rounding of both endpoints unless the enclosure is
                   too wide (then just check containment) */
                decfloat_t r;
                decfloat_init(r, bctx);
                GR_MUST_SUCCEED(decball_set_arb(B, a, bctx));
                GR_MUST_SUCCEED(_decmag_get_decfloat(r, &B->rad, bctx));
                GR_MUST_SUCCEED(decfloat_sub_round(lo, &B->mid, r, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, bctx));
                GR_MUST_SUCCEED(decfloat_add_round(hi, &B->mid, r, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, bctx));
                GR_MUST_SUCCEED(decfloat_set_round(lo, lo, prec, DECIMAL_CTX_RND(ctx), ctx));
                GR_MUST_SUCCEED(decfloat_set_round(hi, hi, prec, DECIMAL_CTX_RND(ctx), ctx));

                if (decfloat_equal(lo, hi, ctx) == T_TRUE && decfloat_equal(lo, y, ctx) != T_TRUE)
                {
                    flint_printf("FAIL: correctly rounded function, which = %d\n", which);
                    gr_ctx_println(ctx);
                    flint_printf("x = %{gr}\n", x, ctx);
                    flint_printf("y = %{gr}\n", y, ctx);
                    flint_printf("lo = %{gr}\n", lo, ctx);
                    flint_printf("a = %{arb}\n", a);
                    flint_abort();
                }

                decfloat_clear(r, bctx);
            }
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(lo, ctx);
        decfloat_clear(hi, ctx);
        decball_clear(B, bctx);
        arb_clear(a);
        gr_ctx_clear(ctx);
        gr_ctx_clear(bctx);
    }

    /* small arguments: the result must equal the rounding of a tight
       high-precision arb enclosure (which needs ~3k digits for x ~ 10^-k) */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, bctx;
        decfloat_t x, y, lo, hi, r;
        decball_t B;
        arb_t a;
        fmpz_t m;
        int which, status, rnd;
        slong prec, k, wpbits;

        gr_ctx_init_decfloat_randtest(ctx, state, 20);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        if (DECIMAL_CTX_IS_EXACT(ctx))
            decimal_ctx_set_prec(ctx, 1 + n_randint(state, 20));
        prec = DECIMAL_CTX_PREC(ctx);
        rnd = DECIMAL_CTX_RND(ctx);

        _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec + 40, DECIMAL_RND_NEAR, 0);

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(lo, ctx);
        decfloat_init(hi, ctx);
        decfloat_init(r, bctx);
        decball_init(B, bctx);
        arb_init(a);
        fmpz_init(m);

        /* x = +/- m 10^-k with 1 <= m < 10^8 and k around the threshold of the shortcut */
        k = prec + 60 + n_randint(state, 200);
        fmpz_set_ui(m, 1 + n_randint(state, 100000000));
        if (n_randint(state, 2)) fmpz_neg(m, m);
        if (n_randint(state, 4) == 0) fmpz_set_si(m, n_randint(state, 2) ? 1 : -1);
        {
            fmpz_t t;
            fmpz_init_set_si(t, -k);
            GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(x, m, t, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
            fmpz_clear(t);
        }

        which = n_randint(state, 11);

        /* log takes 1 + x */
        if (which == 2)
        {
            decfloat_t one;
            decfloat_init(one, ctx);
            GR_MUST_SUCCEED(decfloat_one(one, ctx));
            GR_MUST_SUCCEED(decfloat_add_round(x, x, one, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
            decfloat_clear(one, ctx);
        }

        switch (which)
        {
            case 0: status = decfloat_exp(y, x, ctx); break;
            case 1: status = decfloat_expm1(y, x, ctx); break;
            case 2: status = decfloat_log(y, x, ctx); break;
            case 3: status = decfloat_log1p(y, x, ctx); break;
            case 4: status = decfloat_sin(y, x, ctx); break;
            case 5: status = decfloat_cos(y, x, ctx); break;
            case 6: status = decfloat_tan(y, x, ctx); break;
            case 7: status = decfloat_atan(y, x, ctx); break;
            case 8: status = decfloat_sinh(y, x, ctx); break;
            case 9: status = decfloat_cosh(y, x, ctx); break;
            default: status = decfloat_tanh(y, x, ctx); break;
        }

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: small argument, status %d, which = %d\n", status, which);
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}\n", x, ctx);
            flint_abort();
        }

        wpbits = (slong) (3.33 * (4 * k + prec + 100)) + 64;
        GR_MUST_SUCCEED(decfloat_set_round(&B->mid, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, bctx));
        _decmag_zero(&B->rad, bctx);
        GR_MUST_SUCCEED(decball_get_arb(a, B, wpbits, bctx));

        switch (which)
        {
            case 0: arb_exp(a, a, wpbits); break;
            case 1: arb_expm1(a, a, wpbits); break;
            case 2: arb_log(a, a, wpbits); break;
            case 3: arb_log1p(a, a, wpbits); break;
            case 4: arb_sin(a, a, wpbits); break;
            case 5: arb_cos(a, a, wpbits); break;
            case 6: arb_tan(a, a, wpbits); break;
            case 7: arb_atan(a, a, wpbits); break;
            case 8: arb_sinh(a, a, wpbits); break;
            case 9: arb_cosh(a, a, wpbits); break;
            default: arb_tanh(a, a, wpbits); break;
        }

        decimal_ctx_set_prec(bctx, 4 * k + prec + 100);
        GR_MUST_SUCCEED(decball_set_arb(B, a, bctx));
        GR_MUST_SUCCEED(_decmag_get_decfloat(r, &B->rad, bctx));
        GR_MUST_SUCCEED(decfloat_sub_round(lo, &B->mid, r, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, bctx));
        GR_MUST_SUCCEED(decfloat_add_round(hi, &B->mid, r, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, bctx));
        GR_MUST_SUCCEED(decfloat_set_round(lo, lo, prec, rnd, ctx));
        GR_MUST_SUCCEED(decfloat_set_round(hi, hi, prec, rnd, ctx));

        if (decfloat_equal(lo, hi, ctx) != T_TRUE || decfloat_equal(lo, y, ctx) != T_TRUE)
        {
            flint_printf("FAIL: small argument, which = %d\n", which);
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}\n", x, ctx);
            flint_printf("y = %{gr}\n", y, ctx);
            flint_printf("lo = %{gr}\n", lo, ctx);
            flint_printf("hi = %{gr}\n", hi, ctx);
            flint_abort();
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(lo, ctx);
        decfloat_clear(hi, ctx);
        decfloat_clear(r, bctx);
        decball_clear(B, bctx);
        arb_clear(a);
        fmpz_clear(m);
        gr_ctx_clear(ctx);
        gr_ctx_clear(bctx);
    }

    /* astronomically small arguments where only the shortcut can work:
       the result is x itself or its neighbour, depending on the rounding mode */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, y, ulp, expected;
        fmpz_t t;
        slong prec;
        int rnd, which, status, tail_sign, sgn;

        gr_ctx_init_decfloat_randtest(ctx, state, 20);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        if (DECIMAL_CTX_IS_EXACT(ctx))
            decimal_ctx_set_prec(ctx, 1 + n_randint(state, 20));
        prec = DECIMAL_CTX_PREC(ctx);
        rnd = DECIMAL_CTX_RND(ctx);

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(ulp, ctx);
        decfloat_init(expected, ctx);
        fmpz_init(t);

        /* x = +/- 10^-k with k huge, and possibly more digits */
        fmpz_randtest_unsigned(t, state, 200);
        fmpz_add_ui(t, t, prec + 200);
        fmpz_neg(t, t);
        GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx));
        if (DECFLOAT_IS_SPECIAL(x) || decfloat_digits(x, ctx) > prec)
            GR_MUST_SUCCEED(decfloat_set_si(x, n_randint(state, 2) ? 1 : -1, ctx));
        /* rescale so that the leading digit is at 10^t */
        {
            fmpz_t E;
            fmpz_init(E);
            decfloat_get_sci_exp(E, x, ctx);
            fmpz_sub(t, t, E);
            GR_MUST_SUCCEED(decfloat_mul_10exp_fmpz(x, x, t, ctx));
            fmpz_clear(E);
        }
        sgn = _decfloat_sgn(x, ctx);

        which = n_randint(state, 6);
        switch (which)
        {
            case 0: status = decfloat_sin(y, x, ctx); tail_sign = -sgn; break;
            case 1: status = decfloat_tan(y, x, ctx); tail_sign = sgn; break;
            case 2: status = decfloat_atan(y, x, ctx); tail_sign = -sgn; break;
            case 3: status = decfloat_sinh(y, x, ctx); tail_sign = sgn; break;
            case 4: status = decfloat_tanh(y, x, ctx); tail_sign = -sgn; break;
            default: status = decfloat_expm1(y, x, ctx); tail_sign = 1; break;
        }

        /* expected = round(x + tail_sign * tiny) computed with an explicit tiny */
        decfloat_get_sci_exp(t, x, ctx);
        fmpz_sub_ui(t, t, prec + 3 * DECIMAL_CTX_E(ctx) + 50);
        GR_MUST_SUCCEED(decfloat_set_si(ulp, tail_sign, ctx));
        GR_MUST_SUCCEED(decfloat_mul_10exp_fmpz(ulp, ulp, t, ctx));
        GR_MUST_SUCCEED(decfloat_add_round(expected, x, ulp, prec, rnd, ctx));

        if (status != GR_SUCCESS || decfloat_equal(y, expected, ctx) != T_TRUE)
        {
            flint_printf("FAIL: astronomically small argument, which = %d, status = %d\n", which, status);
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}\n", x, ctx);
            flint_printf("y = %{gr}\n", y, ctx);
            flint_printf("expected = %{gr}\n", expected, ctx);
            flint_abort();
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(ulp, ctx);
        decfloat_clear(expected, ctx);
        fmpz_clear(t);
        gr_ctx_clear(ctx);
    }

    /* exact powers */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, xctx;
        decfloat_t z, x, y, r, expected;
        fmpq_t q;
        slong prec, p, qq;
        int status;

        gr_ctx_init_decfloat_randtest(ctx, state, 20);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        if (DECIMAL_CTX_IS_EXACT(ctx))
            decimal_ctx_set_prec(ctx, 1 + n_randint(state, 20));
        prec = DECIMAL_CTX_PREC(ctx);
        _gr_ctx_init_decimal(xctx, DECIMAL_CTX_FLOAT, DECIMAL_CTX_E(ctx), DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, 0);

        decfloat_init(z, ctx);
        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(r, ctx);
        decfloat_init(expected, ctx);
        fmpq_init(q);

        /* z > 0 with few digits, x = z^qq, y = p/qq with qq | 10^k */
        do { GR_MUST_SUCCEED(decfloat_randtest(z, state, ctx)); } while (DECFLOAT_IS_SPECIAL(z) || decfloat_digits(z, ctx) > 4 || fmpz_bits(&z->exp) > 4);
        z->m.size = FLINT_ABS(z->m.size);
        {
            slong choices[] = {1, 2, 4, 5, 8, 10, 16, 20, 25};
            qq = choices[n_randint(state, 9)];
        }
        p = (slong) n_randint(state, 7) - 3;
        if (p == 0) p = 1;

        GR_MUST_SUCCEED(decfloat_pow_ui(x, z, qq, xctx));
        fmpq_set_si(q, p, qq);
        GR_MUST_SUCCEED(decfloat_set_round_fmpq(y, q, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, xctx));

        status = decfloat_pow(r, x, y, ctx);

        /* expected: z^p rounded */
        if (p > 0)
            GR_MUST_SUCCEED(decfloat_pow_ui(expected, z, p, xctx));
        else
        {
            GR_MUST_SUCCEED(decfloat_pow_ui(expected, z, -p, xctx));
            GR_MUST_SUCCEED(decfloat_get_fmpq(q, expected, xctx));
            fmpq_inv(q, q);
        }
        if (p > 0)
            GR_MUST_SUCCEED(decfloat_set_round(expected, expected, prec, DECIMAL_CTX_RND(ctx), ctx));
        else
            GR_MUST_SUCCEED(decfloat_set_round_fmpq_reference(expected, q, prec, DECIMAL_CTX_RND(ctx), ctx));

        if (status != GR_SUCCESS || decfloat_equal(r, expected, ctx) != T_TRUE)
        {
            flint_printf("FAIL: exact pow, status = %d\n", status);
            gr_ctx_println(ctx);
            flint_printf("z = %{gr}, x = %{gr}, y = %{gr}\n", z, xctx, x, xctx, y, xctx);
            flint_printf("r = %{gr}\n", r, ctx);
            flint_printf("expected = %{gr}\n", expected, ctx);
            flint_abort();
        }

        decfloat_clear(z, ctx);
        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(r, ctx);
        decfloat_clear(expected, ctx);
        fmpq_clear(q);
        gr_ctx_clear(ctx);
        gr_ctx_clear(xctx);
    }

    /* integer powers: correctly rounded, compared with exact rational arithmetic */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, r, expected;
        fmpq_t q;
        fmpz_t n;
        slong prec;
        int status;

        gr_ctx_init_decfloat_randtest(ctx, state, 20);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        if (DECIMAL_CTX_IS_EXACT(ctx))
            decimal_ctx_set_prec(ctx, 1 + n_randint(state, 20));
        prec = DECIMAL_CTX_PREC(ctx);

        decfloat_init(x, ctx);
        decfloat_init(r, ctx);
        decfloat_init(expected, ctx);
        fmpq_init(q);
        fmpz_init(n);

        do { GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx)); } while (DECFLOAT_IS_SPECIAL(x) || decfloat_digits(x, ctx) > 6 || fmpz_bits(&x->exp) > 3);
        fmpz_set_si(n, (slong) n_randint(state, 41) - 20);
        if (n_randint(state, 10) == 0)
            fmpz_set_si(n, (slong) n_randint(state, 2001) - 1000);

        status = decfloat_pow_fmpz(r, x, n, ctx);

        GR_MUST_SUCCEED(decfloat_get_fmpq(q, x, ctx));
        if (fmpz_sgn(n) >= 0)
            fmpq_pow_si(q, q, fmpz_get_si(n));
        else
        {
            fmpq_pow_si(q, q, -fmpz_get_si(n));
            fmpq_inv(q, q);
        }
        GR_MUST_SUCCEED(decfloat_set_round_fmpq_reference(expected, q, prec, DECIMAL_CTX_RND(ctx), ctx));

        if (status != GR_SUCCESS || decfloat_equal(r, expected, ctx) != T_TRUE)
        {
            flint_printf("FAIL: pow_fmpz, status = %d\n", status);
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}, n = %{fmpz}\n", x, ctx, n);
            flint_printf("r = %{gr}\n", r, ctx);
            flint_printf("expected = %{gr}\n", expected, ctx);
            flint_abort();
        }

        decfloat_clear(x, ctx);
        decfloat_clear(r, ctx);
        decfloat_clear(expected, ctx);
        fmpq_clear(q);
        fmpz_clear(n);
        gr_ctx_clear(ctx);
    }

    /* all unary functions: correctly rounded (checked against a tight
       high-precision arb enclosure) for moderate, tiny and special arguments */
    {
        static const struct { int method; int near_one; } tab[] = {
            {GR_METHOD_EXP, 0}, {GR_METHOD_EXPM1, 0}, {GR_METHOD_EXP2, 0},
            {GR_METHOD_EXP10, 0}, {GR_METHOD_LOG, 1}, {GR_METHOD_LOG1P, 0},
            {GR_METHOD_LOG2, 1}, {GR_METHOD_LOG10, 1}, {GR_METHOD_SIN, 0},
            {GR_METHOD_COS, 0}, {GR_METHOD_TAN, 0}, {GR_METHOD_COT, 0},
            {GR_METHOD_SEC, 0}, {GR_METHOD_CSC, 0}, {GR_METHOD_SIN_PI, 0},
            {GR_METHOD_COS_PI, 0}, {GR_METHOD_TAN_PI, 0}, {GR_METHOD_COT_PI, 0},
            {GR_METHOD_CSC_PI, 0}, {GR_METHOD_SEC_PI, 0}, {GR_METHOD_SINC, 0},
            {GR_METHOD_SINC_PI, 0}, {GR_METHOD_ASIN, 0}, {GR_METHOD_ACOS, 0},
            {GR_METHOD_ATAN, 0}, {GR_METHOD_ACOT, 0}, {GR_METHOD_SINH, 0},
            {GR_METHOD_COSH, 0}, {GR_METHOD_TANH, 0}, {GR_METHOD_COTH, 0},
            {GR_METHOD_SECH, 0}, {GR_METHOD_CSCH, 0}, {GR_METHOD_ASINH, 0},
            {GR_METHOD_ACOSH, 1}, {GR_METHOD_ATANH, 0}, {GR_METHOD_ASIN_PI, 0},
            {GR_METHOD_ACOS_PI, 0}, {GR_METHOD_ATAN_PI, 0}, {GR_METHOD_ACOT_PI, 0},
            {GR_METHOD_LAMBERTW, 0}, {GR_METHOD_GAMMA, 0}, {GR_METHOD_RGAMMA, 0},
            {GR_METHOD_LGAMMA, 0}, {GR_METHOD_DIGAMMA, 0}, {GR_METHOD_ZETA, 1},
            {GR_METHOD_ERF, 0}, {GR_METHOD_ERFC, 0}, {GR_METHOD_ERFI, 0},
            {GR_METHOD_ERFINV, 0}, {GR_METHOD_EXP_INTEGRAL_EI, 0}, {GR_METHOD_SIN_INTEGRAL, 0},
            {GR_METHOD_COS_INTEGRAL, 0}, {GR_METHOD_SINH_INTEGRAL, 0}, {GR_METHOD_COSH_INTEGRAL, 0},
            {GR_METHOD_DILOG, 0}, {GR_METHOD_AIRY_AI, 0}, {GR_METHOD_AIRY_BI, 0},
            {GR_METHOD_AGM1, 0}, {GR_METHOD_BARNES_G, 0},
        };
        slong ntab = sizeof(tab) / sizeof(tab[0]);

        for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
        {
            gr_ctx_t ctx, bctx, actx;
            decfloat_t x, y, lo, hi, r;
            decball_t B, X, Y;
            arb_t a, fa;
            fmpz_t m, t;
            slong prec, k, wpbits, which, kind;
            int status, rnd, ok;

            gr_ctx_init_decfloat_randtest(ctx, state, 15);
            decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
            if (DECIMAL_CTX_IS_EXACT(ctx))
                decimal_ctx_set_prec(ctx, 1 + n_randint(state, 15));
            prec = DECIMAL_CTX_PREC(ctx);
            rnd = DECIMAL_CTX_RND(ctx);

            _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec + 40, DECIMAL_RND_NEAR, 0);

            decfloat_init(x, ctx);
            decfloat_init(y, ctx);
            decfloat_init(lo, ctx);
            decfloat_init(hi, ctx);
            decfloat_init(r, bctx);
            decball_init(B, bctx);
            decball_init(X, bctx);
            decball_init(Y, bctx);
            arb_init(a);
            arb_init(fa);
            fmpz_init(m);
            fmpz_init(t);

            which = n_randint(state, ntab);
            kind = n_randint(state, 4);

            switch (kind)
            {
                case 0:
                    /* tiny: +/- m 10^-k */
                    k = prec + 60 + n_randint(state, 100);
                    fmpz_set_ui(m, 1 + n_randint(state, 100000));
                    if (n_randint(state, 4) == 0) fmpz_one(m);
                    if (n_randint(state, 2)) fmpz_neg(m, m);
                    fmpz_set_si(t, -k);
                    GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(x, m, t, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                    break;
                case 1:
                    /* special rational points: k/12, small integers */
                    {
                        fmpq_t q;
                        fmpq_init(q);
                        if (n_randint(state, 2))
                            fmpq_set_si(q, (slong) n_randint(state, 49) - 24, 4);
                        else
                            fmpq_set_si(q, (slong) n_randint(state, 21) - 10, 1 + n_randint(state, 2));
                        GR_MUST_SUCCEED(decfloat_set_round_fmpq(x, q, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                        fmpq_clear(q);
                    }
                    k = 0;
                    break;
                default:
                    /* moderate random */
                    k = 0;
                    GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx));
                    if (!DECFLOAT_IS_SPECIAL(x))
                        fmpz_set_si(&x->exp, (slong) n_randint(state, 4) - 2);
                    if (n_randint(state, 4) == 0)
                        decfloat_neg(x, x, ctx);
                    break;
            }

            if (tab[which].near_one)
            {
                decfloat_t one;
                decfloat_init(one, ctx);
                GR_MUST_SUCCEED(decfloat_one(one, ctx));
                GR_MUST_SUCCEED(decfloat_add_round(x, x, one, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                decfloat_clear(one, ctx);
            }

            status = ((gr_method_unary_op) ctx->methods[tab[which].method])(y, x, ctx);

            wpbits = (slong) (3.33 * (4 * k + 4 * prec + 100)) + 64;
            gr_ctx_init_real_arb(actx, wpbits);

            if (!DECFLOAT_IS_SPECIAL(x))
            {
                GR_MUST_SUCCEED(decfloat_set_round(&B->mid, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, bctx));
                _decmag_zero(&B->rad, bctx);
                GR_MUST_SUCCEED(decball_get_arb(a, B, wpbits, bctx));
            }
            else
                arb_zero(a);

            ok = !DECFLOAT_IS_SPECIAL(x) && (((gr_method_unary_op) actx->methods[tab[which].method])(fa, a, actx) == GR_SUCCESS) && arb_is_finite(fa);

            if (ok)
            {
                decimal_ctx_set_prec(bctx, 4 * k + 4 * prec + 100);
                if (decball_set_arb(Y, fa, bctx) != GR_SUCCESS)
                {
                    flint_printf("FAIL: set_arb of reference, function %wd, kind %wd\n", which, kind);
                    gr_ctx_println(ctx);
                    flint_printf("x = %{gr}\n", x, ctx);
                    flint_printf("a = %{arb}\n", fa);
                    flint_abort();
                }
                GR_MUST_SUCCEED(_decmag_get_decfloat(r, &Y->rad, bctx));
                GR_MUST_SUCCEED(decfloat_sub_round(lo, &Y->mid, r, 4 * k + 4 * prec + 100, DECIMAL_RND_FLOOR, bctx));
                GR_MUST_SUCCEED(decfloat_add_round(hi, &Y->mid, r, 4 * k + 4 * prec + 100, DECIMAL_RND_CEIL, bctx));
                if (decfloat_set_round(lo, lo, prec, rnd, ctx) != GR_SUCCESS || decfloat_set_round(hi, hi, prec, rnd, ctx) != GR_SUCCESS)
                {
                    flint_printf("FAIL: rounding of reference, function %wd, kind %wd\n", which, kind);
                    gr_ctx_println(ctx);
                    flint_printf("x = %{gr}\n", x, ctx);
                    flint_printf("a = %{arb}\n", fa);
                    flint_printf("Y = %{gr}\n", Y, bctx);
                    flint_abort();
                }

                /* the enclosure decides the rounding (which must then agree),
                   or it does not (only possible for an exact value on a
                   boundary, in which case decfloat must still succeed) */
                if (decfloat_equal(lo, hi, ctx) == T_TRUE)
                {
                    if (status != GR_SUCCESS || decfloat_equal(lo, y, ctx) != T_TRUE)
                    {
                        flint_printf("FAIL: unary function %wd, kind %wd, status = %d\n", which, kind, status);
                        gr_ctx_println(ctx);
                        flint_printf("x = %{gr}\n", x, ctx);
                        flint_printf("y = %{gr}\n", y, ctx);
                        flint_printf("expected = %{gr}\n", lo, ctx);
                        flint_printf("a = %{arb}\n", fa);
                        flint_abort();
                    }
                }
                else if (status == GR_SUCCESS && (!_decfloat_is_finite(y) ||
                    _decfloat_cmp(y, lo, ctx) < 0 || _decfloat_cmp(y, hi, ctx) > 0))
                {
                    flint_printf("FAIL: unary function %wd, kind %wd (undecided enclosure), status = %d\n", which, kind, status);
                    gr_ctx_println(ctx);
                    flint_printf("x = %{gr}\n", x, ctx);
                    flint_printf("y = %{gr}\n", y, ctx);
                    flint_printf("lo = %{gr}, hi = %{gr}\n", lo, ctx, hi, ctx);
                    flint_printf("a = %{arb}\n", fa);
                    flint_abort();
                }

                /* the ball function must contain the high-precision value */
                decimal_ctx_set_prec(bctx, prec + 40);
                GR_MUST_SUCCEED(decball_set_decfloat(X, x, bctx));
                if (((gr_method_unary_op) bctx->methods[tab[which].method])(Y, X, bctx) == GR_SUCCESS)
                {
                    decball_t Z;
                    decball_init(Z, bctx);
                    GR_MUST_SUCCEED(decball_set_arb(Z, fa, bctx));
                    if (!_decball_overlaps(Y, Z, bctx))
                    {
                        flint_printf("FAIL: ball function %wd, kind %wd\n", which, kind);
                        gr_ctx_println(bctx);
                        flint_printf("X = %{gr}\n", X, bctx);
                        flint_printf("Y = %{gr}\n", Y, bctx);
                        flint_printf("Z = %{gr}\n", Z, bctx);
                        flint_abort();
                    }
                    decball_clear(Z, bctx);
                }
            }

            gr_ctx_clear(actx);
            decfloat_clear(x, ctx);
            decfloat_clear(y, ctx);
            decfloat_clear(lo, ctx);
            decfloat_clear(hi, ctx);
            decfloat_clear(r, bctx);
            decball_clear(B, bctx);
            decball_clear(X, bctx);
            decball_clear(Y, bctx);
            arb_clear(a);
            arb_clear(fa);
            fmpz_clear(m);
            fmpz_clear(t);
            gr_ctx_clear(ctx);
            gr_ctx_clear(bctx);
        }
    }

    TEST_FUNCTION_END(state);
}
