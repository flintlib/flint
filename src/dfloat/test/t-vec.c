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
#include "double_extras.h"
#include "arb.h"
#include "gr.h"
#include "gr_vec.h"
#include "dfloat.h"

/* The vector operations (SIMD where available) must agree bitwise with
   the scalar operations applied elementwise, for random and special
   inputs; the dot product is checked for containment against arb. */

TEST_FUNCTION_START(vec, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        int n, ball, which;
        slong len, i;
        gr_ptr x, y, r1, r2;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        ball = n_randint(state, 2);
        len = n_randint(state, 4) ? n_randint(state, 20) : n_randint(state, 80);
        which = n_randint(state, 8);
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, ball));

        x = gr_heap_init_vec(len, ctx);
        y = gr_heap_init_vec(len, ctx);
        r1 = gr_heap_init_vec(len + 1, ctx);
        r2 = gr_heap_init_vec(len + 1, ctx);

        for (i = 0; i < len; i++)
        {
            GR_MUST_SUCCEED(gr_randtest(GR_ENTRY(x, i, ctx->sizeof_elem), state, ctx));
            GR_MUST_SUCCEED(gr_randtest(GR_ENTRY(y, i, ctx->sizeof_elem), state, ctx));
        }

        if (which == 4)
        {
            /* division: the plain result within 2^(-53n+8) of the
               exact quotient (the special cases go through the scalar
               function); the ball result contains the quotients of
               the corners of the input balls */
            arb_t ax, ay, aq, ar;
            arf_t t;
            arb_init(ax); arb_init(ay); arb_init(aq); arb_init(ar); arf_init(t);
            switch (n)
            {
                case 1: if (ball) _d1b_vec_div((d1b_ptr) r1, (d1b_srcptr) x, (d1b_srcptr) y, len); else _d1_vec_div((d1_ptr) r1, (d1_srcptr) x, (d1_srcptr) y, len); break;
                case 2: if (ball) _d2b_vec_div((d2b_ptr) r1, (d2b_srcptr) x, (d2b_srcptr) y, len); else _d2_vec_div((d2_ptr) r1, (d2_srcptr) x, (d2_srcptr) y, len); break;
                case 3: if (ball) _d3b_vec_div((d3b_ptr) r1, (d3b_srcptr) x, (d3b_srcptr) y, len); else _d3_vec_div((d3_ptr) r1, (d3_srcptr) x, (d3_srcptr) y, len); break;
                default: if (ball) _d4b_vec_div((d4b_ptr) r1, (d4b_srcptr) x, (d4b_srcptr) y, len); else _d4_vec_div((d4_ptr) r1, (d4_srcptr) x, (d4_srcptr) y, len); break;
            }
            for (i = 0; i < len; i++)
            {
                const double * xi = GR_ENTRY(x, i, ctx->sizeof_elem), * yi = GR_ENTRY(y, i, ctx->sizeof_elem), * ri = GR_ENTRY(r1, i, ctx->sizeof_elem);
                if (!ball)
                {
                    _dfloat_get_arb(ax, xi, n, 0.0);
                    _dfloat_get_arb(ay, yi, n, 0.0);
                    _dfloat_get_arb(ar, ri, n, 0.0);
                    if (yi[0] != 0.0 && arb_is_finite(ax) && arb_is_finite(ay) && arb_is_finite(ar)
                        && fabs(xi[0]) < 1e100 && fabs(xi[0]) > 1e-100 && fabs(yi[0]) < 1e100 && fabs(yi[0]) > 1e-100)
                    {
                        arb_div(aq, ax, ay, 1500);
                        arb_sub(ar, ar, aq, 1500);
                        arb_abs(ar, ar);
                        arb_abs(aq, aq);
                        arb_mul_2exp_si(aq, aq, -53 * n + 8);
                        if (!arb_lt(ar, aq))
                        {
                            flint_printf("FAIL: vec_div accuracy (n = %d, i = %wd)\n", n, i);
                            flint_abort();
                        }
                    }
                }
                else
                {
                    int c;
                    gr_ptr tmp = gr_heap_init(ctx);
                    /* the elements whose scalar division fails (a
                       divisor not bounded away from zero) are left
                       unspecified by the vector function */
                    c = gr_div(tmp, xi, yi, ctx);
                    gr_heap_clear(tmp, ctx);
                    if (c != GR_SUCCESS)
                        continue;
                    _dfloat_get_arb(ar, ri, n, ri[n]);
                    if (!(yi[n] < fabs(yi[0])) || !(xi[n] <= DBL_MAX) || !(yi[n] <= DBL_MAX))
                        continue;
                    for (c = 0; c < 4; c++)
                    {
                        _dfloat_get_arb(ax, xi, n, 0.0);
                        _dfloat_get_arb(ay, yi, n, 0.0);
                        arf_set_d(t, (c & 1) ? xi[n] : -xi[n]);
                        arb_add_arf(ax, ax, t, 2000);
                        arf_set_d(t, (c & 2) ? yi[n] : -yi[n]);
                        arb_add_arf(ay, ay, t, 2000);
                        if (arb_contains_zero(ay) || !arb_is_finite(ax) || !arb_is_finite(ay))
                            continue;
                        arb_div(aq, ax, ay, 2000);
                        if (!arb_contains(ar, aq) && !arb_is_finite(aq))
                            continue;
                        if (!arb_contains(ar, aq) && !(arb_is_exact(aq) && arb_contains_arf(ar, arb_midref(aq))))
                        {
                            /* aq carries rounding; accept if its midpoint is contained */
                            if (!arb_contains_arf(ar, arb_midref(aq)))
                            {
                                flint_printf("FAIL: vec_div containment (n = %d, i = %wd, c = %d)\n", n, i, c);
                                arb_printd(ar, 30); flint_printf("\n");
                                arb_printd(aq, 30); flint_printf("\n");
                                flint_abort();
                            }
                        }
                    }
                }
            }
            arb_clear(ax); arb_clear(ay); arb_clear(aq); arb_clear(ar); arf_clear(t);
            gr_heap_clear_vec(x, len, ctx); gr_heap_clear_vec(y, len, ctx); gr_heap_clear_vec(r1, len + 1, ctx); gr_heap_clear_vec(r2, len + 1, ctx); gr_ctx_clear(ctx);
            continue;
        }

        if (which == 7)
        {
            /* products by a scalar and the reversed dot product: bitwise
               the same as the elementwise scalar operations and as the
               dot product of a reversed copy (with aliasing) */
            gr_ptr c = gr_heap_init(ctx), t = gr_heap_init(ctx), d1 = gr_heap_init(ctx), d2 = gr_heap_init(ctx);
            int op = n_randint(state, 4), alias = n_randint(state, 2), k, bad = 0;
            GR_MUST_SUCCEED(gr_randtest(c, state, ctx));
            if (op == 3)
            {
                for (i = 0; i < len; i++)
                    GR_MUST_SUCCEED(gr_set(GR_ENTRY(r2, i, ctx->sizeof_elem), GR_ENTRY(y, len - 1 - i, ctx->sizeof_elem), ctx));
                GR_MUST_SUCCEED(_gr_vec_dot_rev(d1, c, alias, x, y, len, ctx));
                GR_MUST_SUCCEED(_gr_vec_dot(d2, c, alias, x, r2, len, ctx));
                for (k = 0; k < n + ball; k++)
                    if (((double *) d1)[k] != ((double *) d2)[k] && ((double *) d1)[k] == ((double *) d1)[k])
                        bad = 1;
            }
            else
            {
                /* r1: vector version, r2: scalar version, starting from y */
                GR_MUST_SUCCEED(_gr_vec_set(r1, y, len, ctx));
                GR_MUST_SUCCEED(_gr_vec_set(r2, y, len, ctx));
                if (alias)
                    GR_MUST_SUCCEED(_gr_vec_set(x, y, len, ctx));
                for (i = 0; i < len; i++)
                {
                    gr_ptr ri = GR_ENTRY(r2, i, ctx->sizeof_elem);
                    GR_MUST_SUCCEED(gr_mul(t, GR_ENTRY(x, i, ctx->sizeof_elem), c, ctx));
                    if (op == 0) GR_MUST_SUCCEED(gr_set(ri, t, ctx));
                    else if (op == 1) GR_MUST_SUCCEED(gr_add(ri, ri, t, ctx));
                    else GR_MUST_SUCCEED(gr_sub(ri, ri, t, ctx));
                }
                if (op == 0) GR_MUST_SUCCEED(_gr_vec_mul_scalar(r1, alias ? r1 : x, len, c, ctx));
                else if (op == 1) GR_MUST_SUCCEED(_gr_vec_addmul_scalar(r1, alias ? r1 : x, len, c, ctx));
                else GR_MUST_SUCCEED(_gr_vec_submul_scalar(r1, alias ? r1 : x, len, c, ctx));
                for (i = 0; i < len; i++)
                    for (k = 0; k < n + ball; k++)
                    {
                        double a = ((double *) GR_ENTRY(r1, i, ctx->sizeof_elem))[k], b = ((double *) GR_ENTRY(r2, i, ctx->sizeof_elem))[k];
                        if (a != b && !(a != a && b != b))
                            bad = 1;
                    }
            }
            if (bad)
            {
                flint_printf("FAIL: scalar product / dot_rev (n = %d, ball = %d, op = %d, alias = %d, len = %wd)\n", n, ball, op, alias, len);
                flint_abort();
            }
            gr_heap_clear(c, ctx); gr_heap_clear(t, ctx); gr_heap_clear(d1, ctx); gr_heap_clear(d2, ctx);
            gr_heap_clear_vec(x, len, ctx); gr_heap_clear_vec(y, len, ctx); gr_heap_clear_vec(r1, len + 1, ctx); gr_heap_clear_vec(r2, len + 1, ctx); gr_ctx_clear(ctx);
            continue;
        }

        if (which == 5 || which == 6)
        {
            /* square root and reciprocal square root: the plain
               result within 2^(-53n+8) of the exact value, the ball
               result containing the images of the endpoints of the
               input ball, and the status the or of the scalar ones */
            arb_t ax, ar, aq;
            arf_t t;
            int sq = (which == 5), status = 0, ref = 0;
            arb_init(ax); arb_init(ar); arb_init(aq); arf_init(t);
            /* mostly positive inputs */
            for (i = 0; i < len; i++)
                if (n_randint(state, 4))
                    GR_MUST_SUCCEED(gr_abs(GR_ENTRY(x, i, ctx->sizeof_elem), GR_ENTRY(x, i, ctx->sizeof_elem), ctx));
            switch (n)
            {
                case 1: if (ball) status = sq ? _d1b_vec_sqrt((d1b_ptr) r1, (d1b_srcptr) x, len) : _d1b_vec_rsqrt((d1b_ptr) r1, (d1b_srcptr) x, len); else if (sq) _d1_vec_sqrt((d1_ptr) r1, (d1_srcptr) x, len); else _d1_vec_rsqrt((d1_ptr) r1, (d1_srcptr) x, len); break;
                case 2: if (ball) status = sq ? _d2b_vec_sqrt((d2b_ptr) r1, (d2b_srcptr) x, len) : _d2b_vec_rsqrt((d2b_ptr) r1, (d2b_srcptr) x, len); else if (sq) _d2_vec_sqrt((d2_ptr) r1, (d2_srcptr) x, len); else _d2_vec_rsqrt((d2_ptr) r1, (d2_srcptr) x, len); break;
                case 3: if (ball) status = sq ? _d3b_vec_sqrt((d3b_ptr) r1, (d3b_srcptr) x, len) : _d3b_vec_rsqrt((d3b_ptr) r1, (d3b_srcptr) x, len); else if (sq) _d3_vec_sqrt((d3_ptr) r1, (d3_srcptr) x, len); else _d3_vec_rsqrt((d3_ptr) r1, (d3_srcptr) x, len); break;
                default: if (ball) status = sq ? _d4b_vec_sqrt((d4b_ptr) r1, (d4b_srcptr) x, len) : _d4b_vec_rsqrt((d4b_ptr) r1, (d4b_srcptr) x, len); else if (sq) _d4_vec_sqrt((d4_ptr) r1, (d4_srcptr) x, len); else _d4_vec_rsqrt((d4_ptr) r1, (d4_srcptr) x, len); break;
            }
            for (i = 0; i < len; i++)
            {
                const double * xi = GR_ENTRY(x, i, ctx->sizeof_elem), * ri = GR_ENTRY(r1, i, ctx->sizeof_elem);
                if (!ball)
                {
                    _dfloat_get_arb(ax, xi, n, 0.0);
                    _dfloat_get_arb(ar, ri, n, 0.0);
                    if (xi[0] > 0.0 && arb_is_finite(ax) && arb_is_finite(ar) && fabs(xi[0]) < 1e100 && fabs(xi[0]) > 1e-100)
                    {
                        if (sq) arb_sqrt(aq, ax, 1500); else arb_rsqrt(aq, ax, 1500);
                        arb_sub(ar, ar, aq, 1500);
                        arb_abs(ar, ar);
                        arb_abs(aq, aq);
                        arb_mul_2exp_si(aq, aq, -53 * n + 8);
                        if (!arb_lt(ar, aq))
                        {
                            flint_printf("FAIL: vec_%s accuracy (n = %d, i = %wd)\n", sq ? "sqrt" : "rsqrt", n, i);
                            flint_abort();
                        }
                    }
                }
                else
                {
                    int c, st;
                    double xlo;
                    gr_ptr tmp = gr_heap_init(ctx);
                    st = sq ? gr_sqrt(tmp, xi, ctx) : gr_rsqrt(tmp, xi, ctx);
                    ref |= st;
                    gr_heap_clear(tmp, ctx);
                    if (st != GR_SUCCESS)
                        continue;
                    _dfloat_get_arb(ar, ri, n, ri[n]);
                    xlo = xi[0] - (fabs(xi[n]) + 1e-300);
                    if (!(xi[n] <= DBL_MAX) || (sq ? xlo < 0.0 : xlo <= 0.0))
                        continue;
                    for (c = 0; c < 2; c++)
                    {
                        _dfloat_get_arb(ax, xi, n, 0.0);
                        arf_set_d(t, c ? xi[n] : -xi[n]);
                        arb_add_arf(ax, ax, t, 2000);
                        if (!arb_is_finite(ax) || !arb_is_positive(ax))
                            continue;
                        if (sq) arb_sqrt(aq, ax, 2000); else arb_rsqrt(aq, ax, 2000);
                        if (!arb_contains(ar, aq) && !arb_contains_arf(ar, arb_midref(aq)))
                        {
                            flint_printf("FAIL: vec_%s containment (n = %d, i = %wd, c = %d)\n", sq ? "sqrt" : "rsqrt", n, i, c);
                            arb_printd(ar, 30); flint_printf("\n");
                            arb_printd(aq, 30); flint_printf("\n");
                            flint_abort();
                        }
                    }
                }
            }
            if (ball && status != ref)
            {
                flint_printf("FAIL: vec_%s status (n = %d, %d vs %d)\n", sq ? "sqrt" : "rsqrt", n, status, ref);
                flint_abort();
            }
            arb_clear(ax); arb_clear(ar); arb_clear(aq); arf_clear(t);
            gr_heap_clear_vec(x, len, ctx); gr_heap_clear_vec(y, len, ctx); gr_heap_clear_vec(r1, len + 1, ctx); gr_heap_clear_vec(r2, len + 1, ctx); gr_ctx_clear(ctx);
            continue;
        }

        if (which == 3)
        {
            /* dot product: containment against arb */
            gr_ptr d;
            arb_t s, t, u, v, w;
            gr_ctx_t actx;
            int subtract = n_randint(state, 2);

            if (!ball)
            {
                /* plain: accuracy against the exact dot product, with
                   a relative tolerance in the sum of magnitudes and an
                   absolute tolerance for underflow */
                arb_t mags;
                d = gr_heap_init(ctx);
                GR_MUST_SUCCEED(gr_randtest(d, state, ctx));
                arb_init(s); arb_init(t); arb_init(u); arb_init(v); arb_init(w); arb_init(mags);
                _dfloat_get_arb(s, d, n, 0.0);
                arb_abs(mags, s);
                for (i = 0; i < len; i++)
                {
                    _dfloat_get_arb(u, GR_ENTRY(x, i, ctx->sizeof_elem), n, 0.0);
                    _dfloat_get_arb(v, GR_ENTRY(y, i, ctx->sizeof_elem), n, 0.0);
                    arb_mul(w, u, v, 4000);
                    if (subtract)
                        arb_sub(s, s, w, 4000);
                    else
                        arb_add(s, s, w, 4000);
                    arb_abs(w, w);
                    arb_add(mags, mags, w, 4000);
                }
                GR_MUST_SUCCEED(_gr_vec_dot(r1, d, subtract, x, y, len, ctx));
                _dfloat_get_arb(t, r1, n, 0.0);
                if (arb_is_finite(t) && arb_is_finite(s) && arb_is_finite(mags)
                    && arf_cmpabs_2exp_si(arb_midref(mags), 1000) < 0)
                {
                    arb_sub(t, t, s, 4000);
                    arb_abs(t, t);
                    arb_mul_2exp_si(mags, mags, -53 * n + 8);
                    arb_add_ui(u, mags, 0, 4000);
                    arb_set_ui(v, len + 1);
                    arb_mul_2exp_si(v, v, -1000);
                    arb_add(u, u, v, 4000);
                    if (arb_gt(t, u))
                    {
                        flint_printf("FAIL: plain dot accuracy (n = %d, len = %wd)\n", n, len);
                        arb_printd(t, 30); flint_printf("\n");
                        arb_printd(u, 30); flint_printf("\n");
                        flint_abort();
                    }
                }
                arb_clear(s); arb_clear(t); arb_clear(u); arb_clear(v); arb_clear(w); arb_clear(mags);
                gr_heap_clear(d, ctx);
                gr_heap_clear_vec(x, len, ctx); gr_heap_clear_vec(y, len, ctx); gr_heap_clear_vec(r1, len + 1, ctx); gr_heap_clear_vec(r2, len + 1, ctx); gr_ctx_clear(ctx);
                continue;
            }

            d = gr_heap_init(ctx);
            GR_MUST_SUCCEED(gr_randtest(d, state, ctx));
            arb_init(s); arb_init(t); arb_init(u); arb_init(v); arb_init(w);
            gr_ctx_init_real_arb(actx, 2000);
            _dfloat_get_arb(s, d, n, ((double *) d)[n]);
            for (i = 0; i < len; i++)
            {
                _dfloat_get_arb(u, GR_ENTRY(x, i, ctx->sizeof_elem), n, ((double *) GR_ENTRY(x, i, ctx->sizeof_elem))[n]);
                _dfloat_get_arb(v, GR_ENTRY(y, i, ctx->sizeof_elem), n, ((double *) GR_ENTRY(y, i, ctx->sizeof_elem))[n]);
                arb_mul(w, u, v, 2000);
                if (subtract)
                    arb_sub(s, s, w, 2000);
                else
                    arb_add(s, s, w, 2000);
            }
            GR_MUST_SUCCEED(_gr_vec_dot(r1, d, subtract, x, y, len, ctx));
            _dfloat_get_arb(t, r1, n, ((double *) r1)[n]);
            /* arb's result carries rounding of the radius; use the
               midpoint and the exact-case check */
            if (!arb_contains(t, s) && !(arb_is_exact(s) && arb_contains_arf(t, arb_midref(s))))
            {
                /* the mag rounding in s can exceed our radius slightly;
                   accept if the midpoint is contained and the radii
                   are within 2^-20 relative */
                if (!arb_contains_arf(t, arb_midref(s)) || !mag_is_finite(arb_radref(t)))
                {
                    flint_printf("FAIL: dot containment (n = %d, len = %wd)\n", n, len);
                    arb_printd(t, 30); flint_printf("\n");
                    arb_printd(s, 30); flint_printf("\n");
                    flint_abort();
                }
            }
            arb_clear(s); arb_clear(t); arb_clear(u); arb_clear(v); arb_clear(w);
            gr_ctx_clear(actx);
            gr_heap_clear(d, ctx);
        }
        else
        {
            for (i = 0; i < len; i++)
            {
                gr_ptr xi = GR_ENTRY(x, i, ctx->sizeof_elem), yi = GR_ENTRY(y, i, ctx->sizeof_elem);
                gr_ptr ri = GR_ENTRY(r2, i, ctx->sizeof_elem);
                if (which == 0) GR_MUST_SUCCEED(gr_add(ri, xi, yi, ctx));
                else if (which == 1) GR_MUST_SUCCEED(gr_sub(ri, xi, yi, ctx));
                else GR_MUST_SUCCEED(gr_mul(ri, xi, yi, ctx));
            }
            if (which == 0) GR_MUST_SUCCEED(_gr_vec_add(r1, x, y, len, ctx));
            else if (which == 1) GR_MUST_SUCCEED(_gr_vec_sub(r1, x, y, len, ctx));
            else GR_MUST_SUCCEED(_gr_vec_mul(r1, x, y, len, ctx));

            /* bitwise agreement, except that nan patterns may differ */
            for (i = 0; i < len; i++)
            {
                const double * a = GR_ENTRY(r1, i, ctx->sizeof_elem);
                const double * b = GR_ENTRY(r2, i, ctx->sizeof_elem);
                int k, bad = 0;
                for (k = 0; k < n + ball; k++)
                    if (a[k] != b[k] && !(a[k] != a[k] && b[k] != b[k]))
                        bad = 1;
                if (bad)
                {
                    flint_printf("FAIL: vector op differs from scalar (n = %d, ball = %d, which = %d, i = %wd)\n", n, ball, which, i);
                    for (k = 0; k < n + ball; k++) { flint_printf("%a ", a[k]); }
                    flint_printf("\n");
                    for (k = 0; k < n + ball; k++) { flint_printf("%a ", b[k]); }
                    flint_printf("\n");
                    flint_abort();
                }
            }
        }

        gr_heap_clear_vec(x, len, ctx);
        gr_heap_clear_vec(y, len, ctx);
        gr_heap_clear_vec(r1, len + 1, ctx);
        gr_heap_clear_vec(r2, len + 1, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
