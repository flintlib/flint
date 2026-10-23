/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "double_extras.h"
#include "decimal.h"
#include "fmpq.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"

TEST_FUNCTION_START(decfloat_conv, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, y;
        fmpz_t a, b, m, e;
        fmpq_t q;
        arf_t f, g;
        int status;

        gr_ctx_init_decfloat_randtest(ctx, state, 60);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(m);
        fmpz_init(e);
        fmpq_init(q);
        arf_init(f);
        arf_init(g);

        /* integers roundtrip exactly at exact precision */
        fmpz_randtest(a, state, 1 + n_randint(state, 500));
        status = decfloat_set_round_fmpz(x, a, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);
        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: set_round_fmpz status %d\n", status);
            gr_ctx_println(ctx);
            flint_printf("a = %{fmpz}\n", a);
            flint_abort();
        }
        status = decfloat_get_fmpz(b, x, ctx);
        if (status != GR_SUCCESS || !fmpz_equal(a, b) || !_decfloat_is_int(x, ctx))
        {
            flint_printf("FAIL: fmpz roundtrip\n");
            flint_printf("a = %{fmpz}, x = %{gr}, b = %{fmpz}\n", a, x, ctx, b);
            flint_abort();
        }

        /* m * 10^e decomposition */
        GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx));
        if (DECFLOAT_IS_FINITE(x))
        {
            status = decfloat_get_fmpz_10exp_fmpz(m, e, x, ctx);
            if (status != GR_SUCCESS) { flint_printf("FAIL: get_fmpz_10exp_fmpz %d, x = %{gr}\n", status, x, ctx); gr_ctx_println(ctx); flint_abort(); }
            status = decfloat_set_round_fmpz_10exp_fmpz(y, m, e, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);
            if (status != GR_SUCCESS) { flint_printf("FAIL: set_round_fmpz_10exp_fmpz %d, m = %{fmpz}, e = %{fmpz}\n", status, m, e); gr_ctx_println(ctx); flint_abort(); }
            if (decfloat_equal(x, y, ctx) != T_TRUE || (!fmpz_is_zero(m) && fmpz_divisible_si(m, 10)))
            {
                flint_printf("FAIL: 10exp roundtrip\n");
                flint_printf("x = %{gr}, m = %{fmpz}, e = %{fmpz}, y = %{gr}\n", x, ctx, m, e, y, ctx);
                flint_abort();
            }

            /* fmpq conversion agrees */
            if (decfloat_get_fmpq(q, x, ctx) == GR_SUCCESS)
            {
                GR_MUST_SUCCEED(decfloat_set_round_fmpq(y, q, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                if (decfloat_equal(x, y, ctx) != T_TRUE)
                {
                    flint_printf("FAIL: fmpq roundtrip\n");
                    flint_printf("x = %{gr}, q = %{fmpq}, y = %{gr}\n", x, ctx, q, y, ctx);
                    flint_abort();
                }
            }
        }

        /* arf -> decimal (exact) -> arf */
        arf_randtest(f, state, 1 + n_randint(state, 300), 1 + n_randint(state, 12));
        status = decfloat_set_round_arf(x, f, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);
        if (status == GR_SUCCESS)
        {
            status = decfloat_get_arf(g, x, ARF_PREC_EXACT, ARF_RND_DOWN, ctx);
            if (status != GR_SUCCESS || !arf_equal(f, g))
            {
                flint_printf("FAIL: arf roundtrip\n");
                flint_printf("f = %{arf}, x = %{gr}, g = %{arf}\n", f, x, ctx, g);
                flint_abort();
            }

            /* rounded conversion to arf agrees with arf_set_round of the exact value */
            {
                slong pb = 2 + n_randint(state, 200);
                int rr = n_randint(state, 5);
                arf_t h;
                arf_init(h);
                status = decfloat_get_arf(h, x, pb, rr, ctx); if (status != GR_SUCCESS) { flint_printf("FAIL: get_arf rounded status %d\nf = %{arf}\nx = %{gr}\n", status, f, x, ctx); gr_ctx_println(ctx); flint_abort(); }
                arf_set_round(g, f, pb, rr);
                if (!arf_equal(h, g))
                {
                    flint_printf("FAIL: rounded arf conversion\n");
                    flint_printf("f = %{arf}, x = %{gr}, h = %{arf}, g = %{arf}\n", f, x, ctx, h, g);
                    flint_abort();
                }
                arf_clear(h);
            }
        }

        /* decimal -> arf with huge exponent (Ziv path) vs arb */
        if (n_randint(state, 4) == 0)
        {
            GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx));
            if (!DECFLOAT_IS_SPECIAL(x))
            {
                slong pb = 2 + n_randint(state, 100);
                arb_t t;
                arb_init(t);
                fmpz_set_si(&x->exp, (slong) n_randint(state, 20000) - 10000);
                status = decfloat_get_arf(g, x, pb, ARF_RND_NEAR, ctx);
                if (status == GR_SUCCESS)
                {
                    decball_t bb;
                    gr_ctx_t bctx;
                    _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), 10, DECIMAL_RND_NEAR, 0);
                    decball_init(bb, bctx);
                    GR_MUST_SUCCEED(decfloat_set_round(&bb->mid, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, bctx));
                    GR_MUST_SUCCEED(decball_get_arb(t, bb, pb + 30, bctx));
                    /* g must be within one ulp of the true value */
                    if (!arb_contains_arf(t, g))
                    {
                        arb_t u;
                        arb_init(u);
                        arb_set_arf(u, g);
                        arb_add_error_2exp_si(u, arf_abs_bound_lt_2exp_si(g) - pb);
                        if (!arb_overlaps(u, t))
                        {
                            flint_printf("FAIL: large exponent arf conversion\n");
                            flint_printf("x = %{gr}, g = %{arf}, t = %{arb}\n", x, ctx, g, t);
                            flint_abort();
                        }
                        arb_clear(u);
                    }
                    decball_clear(bb, bctx);
                    gr_ctx_clear(bctx);
                }
                arb_clear(t);
            }
        }

        /* doubles */
        {
            double d = d_randtest(state) * (n_randint(state, 2) ? 1 : -1);
            double d2;
            slong ee = (slong) n_randint(state, 200) - 100;
            d = ldexp(d, ee);
            GR_MUST_SUCCEED(decfloat_set_round_d(x, d, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
            status = decfloat_get_d(&d2, x, ctx); if (status != GR_SUCCESS) { flint_printf("FAIL: get_d status %d, x = %{gr}\n", status, x, ctx); flint_abort(); }
            if (d != d2)
            {
                flint_printf("FAIL: double roundtrip\n");
                flint_printf("d = %.20g, x = %{gr}, d2 = %.20g\n", d, x, ctx, d2);
                flint_abort();
            }
        }

        /* set_other from fmpq and from another decimal context */
        {
            gr_ctx_t QQ, ctx2;
            gr_ctx_init_fmpq(QQ);
            gr_ctx_init_decfloat_randtest(ctx2, state, 60);
            decimal_ctx_set_exp_limits(ctx2, WORD_MIN, WORD_MAX);

            fmpq_randtest(q, state, 100);
            status = gr_set_other(x, q, QQ, ctx);
            if (status == GR_SUCCESS)
            {
                GR_MUST_SUCCEED(decfloat_set_round_fmpq(y, q, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx));
                if (decfloat_equal(x, y, ctx) != T_TRUE)
                {
                    flint_printf("FAIL: set_other fmpq\n");
                    flint_abort();
                }
            }
            else if (!DECIMAL_CTX_IS_EXACT(ctx))
            {
                flint_printf("FAIL: set_other fmpq status %d\n", status);
                flint_abort();
            }

            GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx));
            if (DECFLOAT_IS_FINITE(x))
            {
                decfloat_t z;
                decfloat_init(z, ctx2);
                decimal_ctx_set_prec(ctx2, DECIMAL_PREC_EXACT);
                GR_MUST_SUCCEED(gr_set_other(z, x, ctx, ctx2));
                decimal_ctx_set_prec(ctx, DECIMAL_PREC_EXACT);
                GR_MUST_SUCCEED(gr_set_other(y, z, ctx2, ctx));
                if (decfloat_equal(x, y, ctx) != T_TRUE)
                {
                    flint_printf("FAIL: set_other decfloat\n");
                    flint_printf("x = %{gr}, z = %{gr}, y = %{gr}\n", x, ctx, z, ctx2, y, ctx);
                    flint_abort();
                }
                decfloat_clear(z, ctx2);
            }

            /* through a ball context and back */
            GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx));
            if (DECFLOAT_IS_FINITE(x))
            {
                gr_ctx_t bctx;
                decball_t z;
                decball_t w;

                gr_ctx_init_decball_randtest(bctx, state, 60);
                decimal_ctx_set_exp_limits(bctx, WORD_MIN, WORD_MAX);
                decball_init(z, bctx);
                decball_init(w, bctx);
                decimal_ctx_set_prec(bctx, DECIMAL_PREC_EXACT);
                decimal_ctx_set_prec(ctx, DECIMAL_PREC_EXACT);
                GR_MUST_SUCCEED(gr_set_other(z, x, ctx, bctx));
                GR_MUST_SUCCEED(gr_set_other(y, z, bctx, ctx));
                if (decfloat_equal(x, y, ctx) != T_TRUE || !_decball_is_exact(z, bctx))
                {
                    flint_printf("FAIL: set_other decfloat <-> decball\n");
                    flint_printf("x = %{gr}, z = %{gr}, y = %{gr}\n", x, ctx, z, bctx, y, ctx);
                    flint_abort();
                }

                /* ball to ball with a different limb size or radius precision */
                {
                    gr_ctx_t bctx2;
                    decball_t z2;

                    gr_ctx_init_decball_randtest(bctx2, state, 60);
                    decimal_ctx_set_exp_limits(bctx2, WORD_MIN, WORD_MAX);
                    decball_init(z2, bctx2);

                    decball_add_error_10exp_si(z, -5 + n_randint(state, 10), bctx);
                    GR_MUST_SUCCEED(gr_set_other(z2, z, bctx, bctx2));
                    GR_MUST_SUCCEED(gr_set_other(w, z2, bctx2, bctx));
                    if (!_decball_contains(w, z, bctx))
                    {
                        flint_printf("FAIL: set_other decball <-> decball\n");
                        flint_printf("z = %{gr}, z2 = %{gr}, w = %{gr}\n", z, bctx, z2, bctx2, w, bctx);
                        flint_abort();
                    }

                    decball_clear(z2, bctx2);
                    gr_ctx_clear(bctx2);
                }

                decball_clear(z, bctx);
                decball_clear(w, bctx);
                gr_ctx_clear(bctx);
            }

            gr_ctx_clear(QQ);
            gr_ctx_clear(ctx2);
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(m);
        fmpz_clear(e);
        fmpq_clear(q);
        arf_clear(f);
        arf_clear(g);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
