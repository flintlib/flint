/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"
#include "qqbar.h"
#include "decimal.h"
#include "gr.h"

/* A random real algebraic number: zero, a short decimal number m 10^e
   (sometimes a rounding boundary), a rational number or an irrational
   number of degree 2. Sets digits to the number of significant digits
   of a decimal number and to WORD_MAX otherwise. */
static void
_randtest_part(qqbar_t x, slong * digits, flint_rand_t state, slong prec)
{
    fmpq_t q;
    fmpz_t t;
    int kind = n_randint(state, 4);

    fmpq_init(q);
    fmpz_init(t);
    *digits = WORD_MAX;

    if (kind == 0)
    {
        qqbar_zero(x);
        *digits = 0;
    }
    else if (kind == 1)
    {
        slong e = (slong) n_randint(state, 21) - 10;
        fmpz_randtest_not_zero(fmpq_numref(q), state, 1 + n_randint(state, (ulong) (3.33 * FLINT_MIN(prec, 40) + 10)));
        if (n_randint(state, 2))
            fmpz_mul_ui(fmpq_numref(q), fmpq_numref(q), 10);
        if (n_randint(state, 2))
            fmpz_add_ui(fmpq_numref(q), fmpq_numref(q), 5);
        while (fmpz_divisible_si(fmpq_numref(q), 10) && !fmpz_is_zero(fmpq_numref(q)))
            fmpz_divexact_ui(fmpq_numref(q), fmpq_numref(q), 10);
        *digits = fmpz_is_zero(fmpq_numref(q)) ? 0 : (slong) fmpz_sizeinbase(fmpq_numref(q), 10);
        if (*digits > 1)
        {
            /* sizeinbase may overestimate by one */
            fmpz_ui_pow_ui(t, 10, *digits - 1);
            if (fmpz_cmpabs(fmpq_numref(q), t) < 0)
                (*digits)--;
        }
        fmpz_ui_pow_ui(t, 10, FLINT_ABS(e));
        if (e >= 0)
            fmpz_mul(fmpq_numref(q), fmpq_numref(q), t);
        else
            fmpz_set(fmpq_denref(q), t);
        fmpq_canonicalise(q);
        qqbar_set_fmpq(x, q);
    }
    else if (kind == 2)
    {
        fmpq_randtest_not_zero(q, state, 1 + n_randint(state, 60));
        qqbar_set_fmpq(x, q);
    }
    else
    {
        qqbar_randtest_real(x, state, 2, 1 + n_randint(state, 20));
    }

    if (qqbar_is_zero(x))
        *digits = 0;

    fmpq_clear(q);
    fmpz_clear(t);
}

/* Expected correctly rounded value of the real algebraic number y;
   returns 0 if this is not decided. */
static int
_expected(decfloat_t res, const qqbar_t y, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpq_t q;
    arb_t v;
    gr_ctx_t bctx;
    decball_t Y;
    int ok;

    if (qqbar_degree(y) == 1)
    {
        fmpq_init(q);
        qqbar_get_fmpq(q, y);
        ok = (decfloat_set_round_fmpq(res, q, prec, rnd, ctx) == GR_SUCCESS);
        fmpq_clear(q);
        return ok;
    }

    if (prec == DECIMAL_PREC_EXACT)
        return 0;

    arb_init(v);
    qqbar_get_arb(v, y, 4 * _decimal_digits_to_bits(prec) + 100);
    _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), 2 * prec + 40, DECIMAL_RND_DOWN, 0);
    decimal_ctx_set_rad_prec(bctx, DECMAG_MAX_PREC);
    decball_init(Y, bctx);
    ok = (decball_set_arb(Y, v, bctx) == GR_SUCCESS) &&
         (_decfloat_round_ball(res, Y, prec, rnd, bctx, ctx) == 1) &&
         (_decfloat_finalize(res, ctx) == GR_SUCCESS);
    decball_clear(Y, bctx);
    gr_ctx_clear(bctx);
    arb_clear(v);
    return ok;
}

TEST_FUNCTION_START(decimal_qqbar, state)
{
    slong iter;

    for (iter = 0; iter < 400 * 0.1 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, bctx;
        qqbar_t x, re, im;
        slong prec, dre, dim;
        deccfloat_t z;
        deccball_t b;
        decfloat_t e;
        arb_t v, w;
        int status, k;

        gr_ctx_init_deccfloat_randtest(ctx, state, 30);
        prec = DECIMAL_CTX_PREC(ctx);

        qqbar_init(x);
        qqbar_init(re);
        qqbar_init(im);
        deccfloat_init(z, ctx);
        decfloat_init(e, ctx);
        arb_init(v);
        arb_init(w);

        _randtest_part(re, &dre, state, prec == DECIMAL_PREC_EXACT ? 20 : prec);
        _randtest_part(im, &dim, state, prec == DECIMAL_PREC_EXACT ? 20 : prec);
        qqbar_set_re_im(x, re, im);

        /* correctly rounded parts */
        status = deccfloat_set_qqbar(z, x, ctx);

        for (k = 0; k < 2; k++)
        {
            int rnd = k ? DECIMAL_CTX_RND_IM(ctx) : DECIMAL_CTX_RND(ctx);

            if (status == GR_SUCCESS && _expected(e, k ? im : re, prec, rnd, ctx))
            {
                if (decfloat_equal(e, k ? &z->im : &z->re, ctx) != T_TRUE)
                {
                    flint_printf("FAIL: deccfloat_set_qqbar\n");
                    gr_ctx_println(ctx);
                    flint_printf("x = "); qqbar_printn(x, 30); flint_printf("\n");
                    { char * es = decfloat_get_str(e, ctx); flint_printf("part = %d, expected = %s\n", k, es); flint_free(es); }
                    flint_printf("z = "); gr_println(z, ctx);
                    flint_abort();
                }
            }
        }

        if (qqbar_is_real(x))
        {
            if (_expected(e, re, prec, DECIMAL_CTX_RND(ctx), ctx) &&
                decfloat_set_qqbar(&z->re, x, ctx) == GR_SUCCESS &&
                decfloat_equal(e, &z->re, ctx) != T_TRUE)
            {
                flint_printf("FAIL: decfloat_set_qqbar\n");
                flint_abort();
            }
        }
        else if (decfloat_set_qqbar(e, x, ctx) != GR_DOMAIN)
        {
            flint_printf("FAIL: decfloat_set_qqbar (domain)\n");
            flint_abort();
        }

        /* balls: exact when representable, otherwise enclosures */
        gr_ctx_init_deccball_randtest(bctx, state, 30);
        prec = DECIMAL_CTX_PREC(bctx);
        deccball_init(b, bctx);

        if (deccball_set_qqbar(b, x, bctx) == GR_SUCCESS)
        {
            for (k = 0; k < 2; k++)
            {
                decball_srcptr bp = k ? &b->im : &b->re;
                slong d = k ? dim : dre;

                qqbar_get_arb(v, k ? im : re, 2 * _decimal_digits_to_bits(prec) + 100);
                GR_MUST_SUCCEED(decball_get_arb(w, bp, 2 * _decimal_digits_to_bits(prec) + 100, bctx));

                if (!arb_contains(w, v) || (d <= prec && !DECIMAL_CTX_HAS_EXP_LIMITS(bctx) && !DECMAG_IS_ZERO(&bp->rad)))
                {
                    flint_printf("FAIL: deccball_set_qqbar\n");
                    gr_ctx_println(bctx);
                    flint_printf("x = "); qqbar_printn(x, 30); flint_printf("\n");
                    flint_printf("part = %d, digits = %wd\n", k, d);
                    flint_printf("b = "); gr_println(b, bctx);
                    flint_abort();
                }
            }
        }

        if (qqbar_is_real(x) && decball_set_qqbar(&b->re, x, bctx) == GR_SUCCESS &&
            dre <= prec && !DECIMAL_CTX_HAS_EXP_LIMITS(bctx) && !DECMAG_IS_ZERO(&b->re.rad))
        {
            flint_printf("FAIL: decball_set_qqbar\n");
            flint_abort();
        }

        deccball_clear(b, bctx);
        gr_ctx_clear(bctx);

        qqbar_clear(x);
        qqbar_clear(re);
        qqbar_clear(im);
        deccfloat_clear(z, ctx);
        decfloat_clear(e, ctx);
        arb_clear(v);
        arb_clear(w);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
