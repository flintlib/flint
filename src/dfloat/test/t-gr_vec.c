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
#include "gr.h"
#include "gr_vec.h"
#include "dfloat.h"

/* The generic vector methods through gr on the dfloat rings (and, for
   comparison, the generic fallbacks on arb): the elementary functions
   agree with the scalar methods (overlapping balls; plain values
   within 2^(-53n+10) relative), the conversions between dfloat rings
   contain the input (balls) or are within the target precision
   (plain), gather and scatter agree with elementwise sets, and the
   interval decomposition inverts gr_set_interval_mid_rad. */

static void
vec_get_arb(arb_t res, gr_srcptr x, slong i, int n, int ball, gr_ctx_t ctx)
{
    const double * d = GR_ENTRY(x, i, ctx->sizeof_elem);
    _dfloat_get_arb(res, d, n, ball ? d[n] : 0.0);
}

TEST_FUNCTION_START(gr_vec, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, ctx2, actx;
        int n, n2, ball, ball2, which, strong;
        slong len, i;
        gr_ptr x, y, z, r1, r2, r3;
        arb_t a, b;
        int status;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        n2 = 1 + n_randint(state, DFLOAT_MAX_N);
        ball = n_randint(state, 2);
        ball2 = n_randint(state, 2);
        strong = n_randint(state, 2) ? DFLOAT_STRONG : 0;
        len = n_randint(state, 12);
        which = n_randint(state, 5);
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, ball | strong));
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx2, n2, ball2 | (n_randint(state, 2) ? DFLOAT_STRONG : 0)));
        arb_init(a);
        arb_init(b);

        x = gr_heap_init_vec(len, ctx);
        y = gr_heap_init_vec(len, ctx);
        z = gr_heap_init_vec(len, ctx);
        r1 = gr_heap_init_vec(len, ctx);
        r2 = gr_heap_init_vec(len, ctx);
        r3 = gr_heap_init_vec(len, ctx2);

        for (i = 0; i < len; i++)
        {
            GR_MUST_SUCCEED(gr_randtest(GR_ENTRY(x, i, ctx->sizeof_elem), state, ctx));
            /* moderate arguments for the elementary functions */
            if (which == 0 && n_randint(state, 2))
            {
                arb_randtest(a, state, 53 * n + 10, 6);
                if (n_randint(state, 4) == 0)
                    arb_abs(a, a);
                if (ball)
                    mag_zero(arb_radref(a));
                gr_ctx_init_real_arb(actx, 53 * n + 10);
                GR_MUST_SUCCEED(gr_set_other(GR_ENTRY(x, i, ctx->sizeof_elem), a, actx, ctx));
                gr_ctx_clear(actx);
            }
        }

        if (which == 0)
        {
            /* the elementary functions */
            int f = n_randint(state, 7);
            gr_method_vec_op vecf = NULL;
            gr_method_unary_op scalf = NULL;
            const char * name;

            switch (f)
            {
                case 0: name = "exp"; break;
                case 1: name = "log"; break;
                case 2: name = "sin"; break;
                case 3: name = "cos"; break;
                case 4: name = "sin_cos"; break;
                case 5: name = "sqrt"; break;
                default: name = "rsqrt"; break;
            }
            if (f == 4)
            {
                status = _gr_vec_sin_cos(r1, r2, x, len, ctx);
                for (i = 0; i < len; i++)
                    status |= gr_sin_cos(GR_ENTRY(y, i, ctx->sizeof_elem), GR_ENTRY(z, i, ctx->sizeof_elem), GR_ENTRY(x, i, ctx->sizeof_elem), ctx);
            }
            else
            {
                switch (f)
                {
                    case 0: vecf = GR_VEC_OP(ctx, VEC_EXP); scalf = GR_UNARY_OP(ctx, EXP); break;
                    case 1: vecf = GR_VEC_OP(ctx, VEC_LOG); scalf = GR_UNARY_OP(ctx, LOG); break;
                    case 2: vecf = GR_VEC_OP(ctx, VEC_SIN); scalf = GR_UNARY_OP(ctx, SIN); break;
                    case 3: vecf = GR_VEC_OP(ctx, VEC_COS); scalf = GR_UNARY_OP(ctx, COS); break;
                    case 5: vecf = GR_VEC_OP(ctx, VEC_SQRT); scalf = GR_UNARY_OP(ctx, SQRT); break;
                    default: vecf = GR_VEC_OP(ctx, VEC_RSQRT); scalf = GR_UNARY_OP(ctx, RSQRT); break;
                }
                status = vecf(r1, x, len, ctx);
                for (i = 0; i < len; i++)
                    status |= scalf(GR_ENTRY(y, i, ctx->sizeof_elem), GR_ENTRY(x, i, ctx->sizeof_elem), ctx);
            }

            if (status == GR_SUCCESS)
            {
                for (i = 0; i < len; i++)
                {
                    int k;
                    for (k = 0; k < 1 + (f == 4); k++)
                    {
                        vec_get_arb(a, k ? r2 : r1, i, n, ball, ctx);
                        vec_get_arb(b, k ? z : y, i, n, ball, ctx);
                        if (!arb_is_finite(a) || !arb_is_finite(b))
                            continue;
                        if (ball)
                        {
                            if (!arb_overlaps(a, b))
                            {
                                flint_printf("FAIL: vec_%s overlap (n = %d, i = %wd)\n", name, n, i);
                                arb_printd(a, 30); flint_printf("\n");
                                arb_printd(b, 30); flint_printf("\n");
                                flint_abort();
                            }
                        }
                        else
                        {
                            /* relative, plus absolute in |x| + 1 for
                               the trigonometric functions */
                            arb_sub(a, a, b, 1000);
                            arb_abs(a, a);
                            arb_abs(b, b);
                            if (f >= 2 && f <= 4)
                            {
                                arb_t c;
                                arb_init(c);
                                vec_get_arb(c, x, i, n, ball, ctx);
                                arb_abs(c, c);
                                arb_add_ui(c, c, 1, 1000);
                                arb_add(b, b, c, 1000);
                                arb_clear(c);
                            }
                            arb_mul_2exp_si(b, b, -53 * n + 10);
                            if (!arb_is_zero(a) && !arb_le(a, b))
                            {
                                flint_printf("FAIL: vec_%s accuracy (n = %d, i = %wd)\n", name, n, i);
                                arb_printd(a, 30); flint_printf("\n");
                                arb_printd(b, 30); flint_printf("\n");
                                flint_abort();
                            }
                        }
                    }
                }
            }
        }
        else if (which == 1)
        {
            /* conversions between dfloat rings */
            status = _gr_vec_set_other(r3, x, ctx, len, ctx2);
            if (status != GR_SUCCESS)
            {
                flint_printf("FAIL: vec_set_other status\n");
                flint_abort();
            }
            for (i = 0; i < len; i++)
            {
                vec_get_arb(a, x, i, n, ball, ctx);
                vec_get_arb(b, r3, i, n2, ball2, ctx2);
                if (!arb_is_finite(a) || !arb_is_finite(b))
                    continue;
                if (ball2)
                {
                    /* the target ball contains the source midpoint
                       (source ball: contains the source ball, up to
                       the rounding of the radii to mag_t here) */
                    if (ball)
                    {
                        mag_t t;
                        mag_init(t);
                        mag_mul_2exp_si(t, arb_radref(a), -24);
                        mag_sub_lower(arb_radref(a), arb_radref(a), t);
                        mag_clear(t);
                    }
                    if (ball ? !arb_contains(b, a) : !arb_contains_arf(b, arb_midref(a)))
                    {
                        flint_printf("FAIL: vec_set_other containment (n = %d -> %d)\n", n, n2);
                        arb_printd(a, 30); flint_printf("\n");
                        arb_printd(b, 30); flint_printf("\n");
                        flint_abort();
                    }
                }
                else
                {
                    /* the plain target within the target precision of
                       the source midpoint */
                    mag_zero(arb_radref(a));
                    arb_sub(b, b, a, 1000);
                    arb_abs(b, b);
                    arb_abs(a, a);
                    arb_mul_2exp_si(a, a, -53 * FLINT_MIN(n, n2) + 4);
                    if (!arb_is_zero(b) && !arb_le(b, a))
                    {
                        flint_printf("FAIL: vec_set_other accuracy (n = %d -> %d)\n", n, n2);
                        arb_printd(a, 30); flint_printf("\n");
                        arb_printd(b, 30); flint_printf("\n");
                        flint_abort();
                    }
                }
            }
        }
        else if (which == 2)
        {
            /* gather and scatter against elementwise sets */
            slong * idx = flint_malloc(sizeof(slong) * (len + 1));
            for (i = 0; i < len; i++)
                idx[i] = n_randint(state, len);
            GR_MUST_SUCCEED(_gr_vec_gather(r1, x, idx, len, ctx));
            for (i = 0; i < len; i++)
                GR_MUST_SUCCEED(gr_set(GR_ENTRY(r2, i, ctx->sizeof_elem), GR_ENTRY(x, idx[i], ctx->sizeof_elem), ctx));
            if (memcmp(r1, r2, len * ctx->sizeof_elem) != 0)
            {
                flint_printf("FAIL: vec_gather\n");
                flint_abort();
            }
            GR_MUST_SUCCEED(_gr_vec_set(r1, y, len, ctx));
            GR_MUST_SUCCEED(_gr_vec_set(r2, y, len, ctx));
            /* a permutation, so that the order of the writes is irrelevant */
            for (i = 0; i < len; i++)
                idx[i] = i;
            for (i = 0; i + 1 < len; i++)
            {
                slong j = i + n_randint(state, len - i);
                FLINT_SWAP(slong, idx[i], idx[j]);
            }
            GR_MUST_SUCCEED(_gr_vec_scatter(r1, idx, x, len, ctx));
            for (i = 0; i < len; i++)
                GR_MUST_SUCCEED(gr_set(GR_ENTRY(r2, idx[i], ctx->sizeof_elem), GR_ENTRY(x, i, ctx->sizeof_elem), ctx));
            if (memcmp(r1, r2, len * ctx->sizeof_elem) != 0)
            {
                flint_printf("FAIL: vec_scatter\n");
                flint_abort();
            }
            flint_free(idx);
        }
        else if (which == 3)
        {
            /* the interval decomposition (GR_UNABLE for the whole line,
               whose radius is not a real number) */
            int any_whole = 0, st;
            for (i = 0; i < len && ball; i++)
                any_whole |= _dfloat_ball_is_whole(GR_ENTRY(x, i, ctx->sizeof_elem), n,
                    ((const double *) GR_ENTRY(x, i, ctx->sizeof_elem))[n]);
            st = _gr_vec_get_interval_mid_rad(r1, r2, x, len, ctx);
            if (st != (any_whole ? GR_UNABLE : GR_SUCCESS))
            {
                flint_printf("FAIL: get_interval_mid_rad status\n");
                flint_abort();
            }
            for (i = 0; i < len; i++)
            {
                gr_srcptr m = GR_ENTRY(r1, i, ctx->sizeof_elem), r = GR_ENTRY(r2, i, ctx->sizeof_elem);
                gr_srcptr xi = GR_ENTRY(x, i, ctx->sizeof_elem);
                gr_ptr yi = GR_ENTRY(y, i, ctx->sizeof_elem);
                if (ball && _dfloat_ball_is_whole(xi, n, ((const double *) xi)[n]))
                    continue;
                if (ball)
                {
                    const double * dm = m, * dr = r, * dx = xi;
                    int k;
                    for (k = 0; k < n; k++)
                        if (dm[k] != dx[k] || (k > 0 && dr[k] != 0.0))
                        {
                            flint_printf("FAIL: get_interval_mid_rad components\n");
                            flint_abort();
                        }
                    if (dm[n] != 0.0 || dr[n] != 0.0 || dr[0] != dx[n])
                    {
                        flint_printf("FAIL: get_interval_mid_rad radius\n");
                        flint_abort();
                    }
                    GR_MUST_SUCCEED(gr_set_interval_mid_rad(yi, m, r, ctx));
                    if (!(((const double *) yi)[n] >= dx[n]) || memcmp(yi, xi, n * sizeof(double)) != 0)
                    {
                        flint_printf("FAIL: set_interval_mid_rad inverse\n");
                        flint_abort();
                    }
                }
                else if (!isnan(((const double *) xi)[0]))
                {
                    if (gr_equal(m, xi, ctx) != T_TRUE || gr_is_zero(r, ctx) != T_TRUE)
                    {
                        flint_printf("FAIL: plain get_interval_mid_rad\n");
                        flint_abort();
                    }
                }
            }
        }
        else
        {
            /* the generic fallbacks on arb (not dfloat): exercised for
               coverage of gr_generic */
            gr_ptr ax, ar, as, ac;
            slong * idx = flint_malloc(sizeof(slong) * (len + 1));
            gr_ctx_init_real_arb(actx, 64 + n_randint(state, 100));
            ax = gr_heap_init_vec(len, actx);
            ar = gr_heap_init_vec(len, actx);
            as = gr_heap_init_vec(len, actx);
            ac = gr_heap_init_vec(len, actx);
            for (i = 0; i < len; i++)
            {
                GR_MUST_SUCCEED(gr_randtest(GR_ENTRY(ax, i, actx->sizeof_elem), state, actx));
                idx[i] = n_randint(state, len);
            }
            status = _gr_vec_exp(ar, ax, len, actx);
            status |= _gr_vec_sin_cos(as, ac, ax, len, actx);
            status |= _gr_vec_rsqrt(as, ax, len, actx);
            status |= _gr_vec_get_interval_mid_rad(as, ac, ax, len, actx);
            status |= _gr_vec_gather(ar, ax, idx, len, actx);
            status |= _gr_vec_scatter(as, idx, ar, len, actx);
            status |= _gr_vec_set_other(r1, ax, actx, len, ctx);
            for (i = 0; i < len; i++)
            {
                vec_get_arb(a, r1, i, n, ball, ctx);
                if (ball && arb_is_finite(a) && !arb_contains(a, GR_ENTRY(ax, i, actx->sizeof_elem)))
                {
                    flint_printf("FAIL: vec_set_other from arb\n");
                    flint_abort();
                }
            }
            gr_heap_clear_vec(ax, len, actx);
            gr_heap_clear_vec(ar, len, actx);
            gr_heap_clear_vec(as, len, actx);
            gr_heap_clear_vec(ac, len, actx);
            gr_ctx_clear(actx);
            flint_free(idx);
        }

        gr_heap_clear_vec(x, len, ctx);
        gr_heap_clear_vec(y, len, ctx);
        gr_heap_clear_vec(z, len, ctx);
        gr_heap_clear_vec(r1, len, ctx);
        gr_heap_clear_vec(r2, len, ctx);
        gr_heap_clear_vec(r3, len, ctx2);
        arb_clear(a);
        arb_clear(b);
        gr_ctx_clear(ctx);
        gr_ctx_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
