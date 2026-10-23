/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "fmpq.h"
#include "gr.h"
#include "dfloat.h"

/* the constants (table lookups) against arb: contained (balls) or
   within a relative 2^(2 - 53 N) (plain) */
static void
_t_gr_constants(gr_ctx_t ctx)
{
    int (*f[3])(gr_ptr, gr_ctx_t) = {gr_pi, gr_euler, gr_catalan};
    void (*g[3])(arb_t, slong) = {arb_const_pi, arb_const_euler, arb_const_catalan};
    gr_ctx_t RR;
    gr_ptr x;
    arb_t a, b;
    int i, n = DFLOAT_CTX_N(ctx);
    gr_ctx_init_real_arb(RR, 300);
    x = gr_heap_init(ctx);
    arb_init(a); arb_init(b);
    for (i = 0; i < 3; i++)
    {
        GR_MUST_SUCCEED(f[i](x, ctx));
        GR_MUST_SUCCEED(gr_set_other(a, x, ctx, RR));
        g[i](b, 300);
        if (DFLOAT_CTX_BALL(ctx))
        {
            if (!arb_contains(a, b) || mag_cmp_2exp_si(arb_radref(a), 2 - 53 * n) > 0)
            {
                flint_printf("FAIL: constant %d, n = %d\n", i, n);
                arb_printd(a, 80); flint_printf("\n");
                flint_abort();
            }
        }
        else
        {
            arb_sub(a, a, b, 300);
            arb_div(a, a, b, 300);
            if (arf_cmpabs_2exp_si(arb_midref(a), 2 - 53 * n) > 0)
            {
                flint_printf("FAIL: constant %d, n = %d\n", i, n);
                flint_abort();
            }
        }
    }
    arb_clear(a); arb_clear(b);
    gr_heap_clear(x, ctx);
    gr_ctx_clear(RR);
}

/* conversions to acb agree with those to arb (the imaginary part
   zero), including the nonfinite values (GR_DOMAIN for the plain types,
   [+/- inf] for the balls) */
static void
_t_gr_to_acb(gr_ctx_t ctx, flint_rand_t state)
{
    gr_ctx_t RR, CC;
    gr_ptr x;
    arb_t a;
    acb_t z;
    int i, st1, st2;

    gr_ctx_init_real_arb(RR, 300);
    gr_ctx_init_complex_acb(CC, 300);
    x = gr_heap_init(ctx);
    arb_init(a);
    acb_init(z);
    for (i = 0; i < 100; i++)
    {
        GR_MUST_SUCCEED(gr_randtest(x, state, ctx));
        st1 = gr_set_other(a, x, ctx, RR);
        st2 = gr_set_other(z, x, ctx, CC);
        if (st1 != st2 || (st1 == GR_SUCCESS &&
            (!arb_equal(acb_realref(z), a) || !arb_is_zero(acb_imagref(z)))))
        {
            flint_printf("FAIL: set_other (acb), n = %d\n", DFLOAT_CTX_N(ctx));
            gr_println(x, ctx);
            arb_printd(a, 30); flint_printf("\n");
            acb_printd(z, 30); flint_printf("\n");
            flint_abort();
        }
    }
    arb_clear(a);
    acb_clear(z);
    gr_heap_clear(x, ctx);
    gr_ctx_clear(RR);
    gr_ctx_clear(CC);
}

/* the integer and rational interface on values that the type
   represents exactly (|v| < 2^50): set, get, cmp, cmpabs, abs, sgn,
   and the strings */
static void
_t_gr_integers(gr_ctx_t ctx, flint_rand_t state)
{
    gr_ptr x, y, z;
    fmpz_t f, g;
    fmpq_t q;
    slong v, w, sv;
    ulong uv;
    int i, c, n = DFLOAT_CTX_N(ctx);
    char * str;

    x = gr_heap_init(ctx);
    y = gr_heap_init(ctx);
    z = gr_heap_init(ctx);
    fmpz_init(f);
    fmpz_init(g);
    fmpq_init(q);
    for (i = 0; i < 100; i++)
    {
#if FLINT_BITS == 64
        v = (slong) n_randint(state, UWORD(1) << 25) * (slong) n_randint(state, UWORD(1) << 24);
#else
        v = (slong) n_randint(state, UWORD(1) << 30);
#endif
        if (n_randint(state, 2))
            v = -v;
        w = (slong) n_randint(state, 1000) - 500;

        GR_MUST_SUCCEED(gr_set_si(x, v, ctx));
        GR_MUST_SUCCEED(gr_set_si(y, w, ctx));
        sv = 0;
        if (gr_get_si(&sv, x, ctx) != GR_SUCCESS || sv != v)
        {
            flint_printf("FAIL: set_si / get_si, n = %d, v = %wd\n", n, v);
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_get_fmpz(f, x, ctx));
        if (!fmpz_equal_si(f, v))
        {
            flint_printf("FAIL: get_fmpz, n = %d\n", n);
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_set_fmpz(z, f, ctx));
        GR_MUST_SUCCEED(gr_abs(y, x, ctx));
        GR_MUST_SUCCEED(gr_get_fmpz(g, y, ctx));
        uv = 0;
        if (gr_equal(z, x, ctx) != T_TRUE || fmpz_cmpabs(g, f) != 0 || fmpz_sgn(g) < 0 ||
            gr_get_ui(&uv, y, ctx) != GR_SUCCESS || uv != fmpz_get_ui(g))
        {
            flint_printf("FAIL: set_fmpz / abs / get_ui, n = %d\n", n);
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_set_ui(z, uv, ctx));
        GR_MUST_SUCCEED(gr_sgn(y, x, ctx));
        GR_MUST_SUCCEED(gr_get_si(&sv, y, ctx));
        if (gr_equal(z, x, ctx) == (v < 0 ? T_TRUE : T_FALSE) || sv != (v > 0) - (v < 0))
        {
            flint_printf("FAIL: set_ui / sgn, n = %d\n", n);
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_set_si(y, w, ctx));
        if (gr_cmp(&c, x, y, ctx) != GR_SUCCESS || c != (v > w) - (v < w) ||
            gr_cmpabs(&c, x, y, ctx) != GR_SUCCESS || c != (FLINT_ABS(v) > FLINT_ABS(w)) - (FLINT_ABS(v) < FLINT_ABS(w)))
        {
            flint_printf("FAIL: cmp / cmpabs, n = %d, v = %wd, w = %wd\n", n, v, w);
            flint_abort();
        }
        /* v / 4 exactly, through fmpq */
        fmpq_set_si(q, v, 4);
        GR_MUST_SUCCEED(gr_set_fmpq(y, q, ctx));
        GR_MUST_SUCCEED(gr_mul_2exp_si(z, x, -2, ctx));
        if (gr_equal(y, z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: set_fmpq, n = %d\n", n);
            flint_abort();
        }
        /* strings: the integer prints and parses back */
        GR_MUST_SUCCEED(gr_get_str(&str, x, ctx));
        if (gr_set_str(y, str, ctx) == GR_SUCCESS && gr_equal(y, x, ctx) == T_FALSE)
        {
            flint_printf("FAIL: get_str / set_str, n = %d, %s\n", n, str);
            flint_abort();
        }
        flint_free(str);
    }
    gr_heap_clear(x, ctx);
    gr_heap_clear(y, ctx);
    gr_heap_clear(z, ctx);
    fmpz_clear(f);
    fmpz_clear(g);
    fmpq_clear(q);
}

TEST_FUNCTION_START(gr, state)
{
    gr_ctx_t ctx;
    int n, flags;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (n = 1; n <= DFLOAT_MAX_N; n++)
    {
        for (flags = 0; flags <= 7; flags++)
        {
            int ball = flags & DFLOAT_BALL;
            GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, flags));
            _t_gr_constants(ctx);
            _t_gr_to_acb(ctx, state);
            _t_gr_integers(ctx, state);
            if (ball)
                gr_test_ring(ctx, 300, 0);
            else
                gr_test_floating_point(ctx, 300, 0);
            gr_ctx_clear(ctx);
        }
    }

    TEST_FUNCTION_END(state);
}
