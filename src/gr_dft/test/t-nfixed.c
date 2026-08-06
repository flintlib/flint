/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "arf.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* random arf in (-1, 1) */
static void
nfixed_randarf(arf_t x, flint_rand_t state, slong nlimbs)
{
    arf_randtest(x, state, nlimbs * FLINT_BITS, 0);
    while (arf_cmpabs_2exp_si(x, 0) >= 0)
        arf_mul_2exp_si(x, x, -1);
}

TEST_FUNCTION_START(gr_dft_nfixed, state)
{
    slong iter;

    /* arithmetic against arf reference */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong nlimbs = 1 + n_randint(state, 24);
        gr_ctx_t ctx;
        gr_ptr a, b, c;
        arf_t fa, fb, fc, ref, ulp;
        int status = GR_SUCCESS;

        if (gr_dft_ctx_init_nfixed(ctx, nlimbs) != GR_SUCCESS)
            flint_abort();

        a = gr_heap_init(ctx);
        b = gr_heap_init(ctx);
        c = gr_heap_init(ctx);
        arf_init(fa); arf_init(fb); arf_init(fc);
        arf_init(ref); arf_init(ulp);

        nfixed_randarf(fa, state, nlimbs);
        nfixed_randarf(fb, state, nlimbs);
        status |= gr_dft_nfixed_set_arf(a, fa, ctx);
        status |= gr_dft_nfixed_set_arf(b, fb, ctx);
        gr_dft_nfixed_get_arf(fa, a, ctx);
        gr_dft_nfixed_get_arf(fb, b, ctx);

        /* multiplication: at most 2 ulp error */
        status |= gr_mul(c, a, b, ctx);
        gr_dft_nfixed_get_arf(fc, c, ctx);
        arf_mul(ref, fa, fb, 2 * nlimbs * FLINT_BITS + 32, ARF_RND_NEAR);
        arf_sub(ref, ref, fc, 2 * nlimbs * FLINT_BITS + 32, ARF_RND_NEAR);
        arf_abs(ref, ref);
        arf_one(ulp);
        arf_mul_2exp_si(ulp, ulp, -nlimbs * FLINT_BITS + 1);
        if (arf_cmp(ref, ulp) > 0)
            TEST_FUNCTION_FAIL("mul error exceeds 2 ulp, nlimbs = %wd\n", nlimbs);

        /* complex multiplication: at most cmul_err ulp per component,
           covering the sized, schoolbook and Karatsuba paths */
        {
            gr_ctx_t cctx;
            gr_ptr ca, cb, cc;
            arf_t gr, gi, rr, ri;
            slong prec = 2 * nlimbs * FLINT_BITS + 32;

            if (gr_dft_ctx_init_nfixed_complex(cctx, nlimbs) != GR_SUCCESS)
                flint_abort();

            ca = gr_heap_init(cctx);
            cb = gr_heap_init(cctx);
            cc = gr_heap_init(cctx);
            arf_init(gr); arf_init(gi); arf_init(rr); arf_init(ri);

            /* local copies of fa, fb halved as the real parts (a, b,
               fa, fb are reused by the checks below and must not be
               clobbered); fresh imaginary parts, all halved so no
               component of the product can reach 1 */
            arf_t ha, hb;
            arf_init(ha);
            arf_init(hb);
            arf_mul_2exp_si(ha, fa, -1);
            arf_mul_2exp_si(hb, fb, -1);
            nfixed_randarf(gr, state, nlimbs);
            arf_mul_2exp_si(gr, gr, -1);
            nfixed_randarf(gi, state, nlimbs);
            arf_mul_2exp_si(gi, gi, -1);

            status |= gr_dft_nfixed_set_arf(ca, ha, ctx);
            status |= gr_dft_nfixed_set_arf(
                    (gr_ptr)((char *) ca + ctx->sizeof_elem), gr, ctx);
            status |= gr_dft_nfixed_set_arf(cb, hb, ctx);
            status |= gr_dft_nfixed_set_arf(
                    (gr_ptr)((char *) cb + ctx->sizeof_elem), gi, ctx);
            gr_dft_nfixed_get_arf(ha, ca, ctx);
            gr_dft_nfixed_get_arf(gr,
                    (gr_srcptr)((const char *) ca + ctx->sizeof_elem), ctx);
            gr_dft_nfixed_get_arf(hb, cb, ctx);
            gr_dft_nfixed_get_arf(gi,
                    (gr_srcptr)((const char *) cb + ctx->sizeof_elem), ctx);

            status |= gr_mul(cc, ca, cb, cctx);

            arf_mul(rr, ha, hb, prec, ARF_RND_NEAR);
            arf_mul(ref, gr, gi, prec, ARF_RND_NEAR);
            arf_sub(rr, rr, ref, prec, ARF_RND_NEAR);
            arf_mul(ri, ha, gi, prec, ARF_RND_NEAR);
            arf_mul(ref, gr, hb, prec, ARF_RND_NEAR);
            arf_add(ri, ri, ref, prec, ARF_RND_NEAR);

            arf_set_d(ulp, _gr_dft_nfixed_cmul_err_ulps(cctx));
            arf_mul_2exp_si(ulp, ulp, -nlimbs * FLINT_BITS);

            gr_dft_nfixed_get_arf(fc, cc, ctx);
            arf_sub(ref, fc, rr, prec, ARF_RND_UP);
            arf_abs(ref, ref);
            if (arf_cmp(ref, ulp) > 0)
                TEST_FUNCTION_FAIL("cmul re error, nlimbs = %wd\n", nlimbs);
            gr_dft_nfixed_get_arf(fc,
                    (gr_srcptr)((const char *) cc + ctx->sizeof_elem), ctx);
            arf_sub(ref, fc, ri, prec, ARF_RND_UP);
            arf_abs(ref, ref);
            if (arf_cmp(ref, ulp) > 0)
                TEST_FUNCTION_FAIL("cmul im error, nlimbs = %wd\n", nlimbs);

            gr_heap_clear(ca, cctx);
            gr_heap_clear(cb, cctx);
            gr_heap_clear(cc, cctx);
            arf_clear(ha); arf_clear(hb);
            arf_clear(gr); arf_clear(gi); arf_clear(rr); arf_clear(ri);
            gr_ctx_clear(cctx);
        }

        /* addition: exact and GR_SUCCESS when the result is in
           range; GR_UNABLE exactly when the exact result has
           magnitude at least 1 (a successful status certifies that
           no overflow occurred) */
        {
            int st = gr_add(c, a, b, ctx);
            arf_add(ref, fa, fb, 2 * nlimbs * FLINT_BITS + 32,
                    ARF_RND_NEAR);
            if (arf_cmpabs_2exp_si(ref, 0) < 0)
            {
                gr_dft_nfixed_get_arf(fc, c, ctx);
                if (st != GR_SUCCESS || !arf_equal(ref, fc))
                    TEST_FUNCTION_FAIL("add not exact, nlimbs = %wd\n",
                            nlimbs);
            }
            else if (st != GR_UNABLE)
                TEST_FUNCTION_FAIL("add overflow not flagged, "
                        "nlimbs = %wd\n", nlimbs);
        }

        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("unexpected status\n");

        gr_heap_clear(a, ctx);
        gr_heap_clear(b, ctx);
        gr_heap_clear(c, ctx);
        arf_clear(fa); arf_clear(fb); arf_clear(fc);
        arf_clear(ref); arf_clear(ulp);
        gr_ctx_clear(ctx);
    }

    /* DFT against a high-precision acb reference of the exact stored
       inputs; the observed errors must respect the bound function */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        static const slong nl_tab[] = { 1, 2, 3, 4, 6, 16 };
        slong nlimbs = nl_tab[n_randint(state, 6)];
        ulong n = 1 + n_randint(state, 80);
        gr_ctx_t rctx, cctx, actx;
        gr_dft_pre_t P, PA;
        gr_ptr x, y, xa, ya;
        double peak1, err1, peak, err, scale;
        slong prec = nlimbs * FLINT_BITS + 128;
        slong asz;
        ulong j;
        arf_t t, s;
        int status = GR_SUCCESS;

        int alg = GR_DFT_ALG_AUTO;

        if (n_randint(state, 2))
            n = UWORD(1) << (1 + n_randint(state, 9));
        else if (n % 2 == 1 && n >= 3 && n_randint(state, 2))
            alg = GR_DFT_ALG_BLUESTEIN;

        if (gr_dft_ctx_init_nfixed(rctx, nlimbs) != GR_SUCCESS)
            flint_abort();
        if (gr_dft_ctx_init_nfixed_complex(cctx, nlimbs) != GR_SUCCESS)
            flint_abort();
        gr_ctx_init_complex_acb(actx, prec);
        asz = actx->sizeof_elem;

        status |= gr_dft_precomp_init_karatsuba(P, n, alg, 0, rctx, cctx);
        status |= gr_dft_precomp_init_karatsuba(PA, n, alg, 0, NULL, actx);
        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("plan init failed, n = %wu, nlimbs = %wd\n",
                    n, nlimbs);

        gr_dft_precomp_nfixed_bound(&peak1, &err1, 1.0, 0.0, P);
        scale = 0.5 / peak1;

        x = gr_heap_init_vec(n, cctx);
        y = gr_heap_init_vec(n, cctx);
        xa = gr_heap_init_vec(n, actx);
        ya = gr_heap_init_vec(n, actx);
        arf_init(t);
        arf_init(s);

        for (j = 0; j < 2 * n; j++)
        {
            gr_ptr comp = (gr_ptr) ((char *) x + j * rctx->sizeof_elem);
            acb_ptr aj = (acb_ptr) GR_ENTRY(xa, j / 2, asz);
            arb_ptr part = (j % 2) ? acb_imagref(aj) : acb_realref(aj);

            nfixed_randarf(t, state, nlimbs);
            arf_set_d(s, scale);
            arf_mul(t, t, s, prec, ARF_RND_DOWN);
            status |= gr_dft_nfixed_set_arf(comp, t, rctx);
            gr_dft_nfixed_get_arf(arb_midref(part), comp, rctx);
        }

        status |= gr_dft_precomp(y, x, P, cctx);
        status |= gr_dft_precomp(ya, xa, PA, actx);
        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("transform failed, n = %wu, nlimbs = %wd\n",
                    n, nlimbs);

        gr_dft_precomp_nfixed_bound(&peak, &err, scale, 0.0, P);
        if (peak >= 1.0)
            TEST_FUNCTION_FAIL("peak >= 1 despite scaling, n = %wu\n", n);

        for (j = 0; j < 2 * n; j++)
        {
            gr_srcptr comp = (gr_srcptr)
                ((const char *) y + j * rctx->sizeof_elem);
            acb_srcptr aj = (acb_srcptr) GR_ENTRY(ya, j / 2, asz);
            arb_srcptr part = (j % 2) ? acb_imagref(aj) : acb_realref(aj);
            arf_t d, bnd;

            arf_init(d);
            arf_init(bnd);

            gr_dft_nfixed_get_arf(t, comp, rctx);
            arf_sub(d, t, arb_midref(part), prec, ARF_RND_UP);
            arf_abs(d, d);
            arf_set_mag(bnd, arb_radref(part));
            arf_add(d, d, bnd, prec, ARF_RND_UP);

            arf_set_d(bnd, err);
            arf_mul_2exp_si(bnd, bnd, -nlimbs * FLINT_BITS);

            if (arf_cmp(d, bnd) > 0)
                TEST_FUNCTION_FAIL("error exceeds bound, n = %wu, "
                        "nlimbs = %wd, component %wu\n", n, nlimbs, j);

            arf_clear(d);
            arf_clear(bnd);
        }

        arf_clear(t);
        arf_clear(s);
        gr_heap_clear_vec(x, n, cctx);
        gr_heap_clear_vec(y, n, cctx);
        gr_heap_clear_vec(xa, n, actx);
        gr_heap_clear_vec(ya, n, actx);
        gr_dft_precomp_clear(P);
        gr_dft_precomp_clear(PA);
        gr_ctx_clear(rctx);
        gr_ctx_clear(cctx);
        gr_ctx_clear(actx);
    }


    /* Method-table surface not reached through the transforms:
       constants, integer setters, squaring, predicates, random
       generation, the exported complex multiplication kernels, and
       context printing. */
    {
        slong sizes[] = { 1, 2, 3, 4, 6, 24 };
        slong si;

        for (si = 0; si < 6; si++)
        {
            slong nlimbs = sizes[si];
            gr_ctx_t ctx, cctx;
            gr_ptr a, b, c, d;
            arf_t fa, fb, ref;
            char * str = NULL;

            GR_MUST_SUCCEED(gr_dft_ctx_init_nfixed(ctx, nlimbs));
            GR_MUST_SUCCEED(gr_dft_ctx_init_nfixed_complex(cctx, nlimbs));
            a = gr_heap_init(ctx);
            b = gr_heap_init(ctx);
            arf_init(fa); arf_init(fb); arf_init(ref);

            /* one and neg_one clamp to 1 - ulp and must cancel */
            GR_MUST_SUCCEED(gr_one(a, ctx));
            GR_MUST_SUCCEED(gr_neg_one(b, ctx));
            gr_dft_nfixed_get_arf(fa, a, ctx);
            gr_dft_nfixed_get_arf(fb, b, ctx);
            arf_add(ref, fa, fb, nlimbs * FLINT_BITS + 8, ARF_RND_NEAR);
            if (!arf_is_zero(ref) || gr_is_zero(a, ctx) != T_FALSE)
                TEST_FUNCTION_FAIL("one/neg_one, nlimbs = %wd\n", nlimbs);

            /* integer setters: 0 and +-1 in range, others GR_DOMAIN */
            if (gr_set_si(a, 0, ctx) != GR_SUCCESS ||
                gr_is_zero(a, ctx) != T_TRUE ||
                gr_set_si(a, 1, ctx) != GR_SUCCESS ||
                gr_set_si(a, -1, ctx) != GR_SUCCESS ||
                gr_set_si(a, 2, ctx) != GR_DOMAIN ||
                gr_set_ui(a, 1, ctx) != GR_SUCCESS ||
                gr_set_ui(a, 7, ctx) != GR_DOMAIN)
                TEST_FUNCTION_FAIL("set_si/set_ui, nlimbs = %wd\n", nlimbs);

            /* squaring: within 2 ulp of the exact square, and equal
               semantics on the operands */
            GR_MUST_SUCCEED(gr_randtest(a, state, ctx));
            GR_MUST_SUCCEED(gr_set(b, a, ctx));
            if (gr_equal(a, b, ctx) != T_TRUE)
                TEST_FUNCTION_FAIL("equal, nlimbs = %wd\n", nlimbs);
            GR_MUST_SUCCEED(gr_sqr(b, a, ctx));
            gr_dft_nfixed_get_arf(fa, a, ctx);
            gr_dft_nfixed_get_arf(fb, b, ctx);
            arf_mul(ref, fa, fa, 2 * nlimbs * FLINT_BITS + 8, ARF_RND_NEAR);
            arf_sub(ref, ref, fb, 2 * nlimbs * FLINT_BITS + 8, ARF_RND_NEAR);
            arf_abs(ref, ref);
            arf_mul_2exp_si(ref, ref, nlimbs * FLINT_BITS - 1);
            if (arf_cmp_si(ref, 1) > 0)
                TEST_FUNCTION_FAIL("sqr error > 2 ulp, nlimbs = %wd\n",
                        nlimbs);

            /* aliased generic multiplication and squaring stage
               through a temporary; swap for completeness */
            GR_MUST_SUCCEED(gr_randtest(a, state, ctx));
            GR_MUST_SUCCEED(gr_randtest(b, state, ctx));
            GR_MUST_SUCCEED(gr_mul(a, a, b, ctx));
            GR_MUST_SUCCEED(gr_sqr(a, a, ctx));
            gr_swap(a, b, ctx);

            /* randtest stays in range */
            GR_MUST_SUCCEED(gr_randtest(a, state, ctx));
            gr_dft_nfixed_get_arf(fa, a, ctx);
            if (arf_cmpabs_2exp_si(fa, 0) >= 0)
                TEST_FUNCTION_FAIL("randtest out of range, nlimbs = %wd\n",
                        nlimbs);

            /* exported complex multiplication kernels agree within
               their documented error bound */
            c = gr_heap_init(cctx);
            d = gr_heap_init(cctx);
            GR_MUST_SUCCEED(gr_randtest(c, state, cctx));
            GR_MUST_SUCCEED(gr_randtest(d, state, cctx));
            {
                gr_ptr r1 = gr_heap_init(cctx), r2 = gr_heap_init(cctx);
                arf_t g1, g2, dd;
                double eb = _gr_dft_nfixed_cmul_err_ulps(cctx);
                slong part;

                GR_MUST_SUCCEED(_gr_dft_nfixed_cmul_schoolbook(r1, c, d, cctx));
                GR_MUST_SUCCEED(_gr_dft_nfixed_cmul_karatsuba(r2, c, d, cctx));

                arf_init(g1); arf_init(g2); arf_init(dd);
                for (part = 0; part < 2; part++)
                {
                    gr_dft_nfixed_get_arf(g1, (gr_srcptr)
                            ((const char *) r1 + part * (nlimbs + 1) *
                             (slong) sizeof(ulong)), ctx);
                    gr_dft_nfixed_get_arf(g2, (gr_srcptr)
                            ((const char *) r2 + part * (nlimbs + 1) *
                             (slong) sizeof(ulong)), ctx);
                    arf_sub(dd, g1, g2, 2 * nlimbs * FLINT_BITS, ARF_RND_NEAR);
                    arf_abs(dd, dd);
                    arf_mul_2exp_si(dd, dd, nlimbs * FLINT_BITS);
                    if (arf_cmp_d(dd, 2 * eb) > 0)
                        TEST_FUNCTION_FAIL("cmul kernels disagree, "
                                "nlimbs = %wd\n", nlimbs);
                }
                arf_clear(g1); arf_clear(g2); arf_clear(dd);
                gr_heap_clear(r1, cctx);
                gr_heap_clear(r2, cctx);
            }
            gr_heap_clear(c, cctx);
            gr_heap_clear(d, cctx);

            /* the same surface on the complex context */
            c = gr_heap_init(cctx);
            GR_MUST_SUCCEED(gr_one(c, cctx));
            GR_MUST_SUCCEED(gr_neg(c, c, cctx));
            d = gr_heap_init(cctx);
            GR_MUST_SUCCEED(gr_neg_one(d, cctx));
            if (gr_equal(c, d, cctx) != T_TRUE ||
                gr_is_zero(c, cctx) != T_FALSE)
                TEST_FUNCTION_FAIL("complex constants, nlimbs = %wd\n",
                        nlimbs);
            if (gr_set_si(c, 0, cctx) != GR_SUCCESS ||
                gr_is_zero(c, cctx) != T_TRUE ||
                gr_set_si(c, -1, cctx) != GR_SUCCESS ||
                gr_set_si(c, 3, cctx) != GR_DOMAIN ||
                gr_set_ui(c, 1, cctx) != GR_SUCCESS ||
                gr_set_ui(c, 9, cctx) != GR_DOMAIN)
                TEST_FUNCTION_FAIL("complex setters, nlimbs = %wd\n",
                        nlimbs);
            GR_MUST_SUCCEED(gr_randtest(c, state, cctx));
            GR_MUST_SUCCEED(gr_set(d, c, cctx));
            if (gr_equal(c, d, cctx) != T_TRUE)
                TEST_FUNCTION_FAIL("complex equal, nlimbs = %wd\n", nlimbs);
            GR_MUST_SUCCEED(gr_sqr(d, c, cctx));
            gr_heap_clear(c, cctx);
            gr_heap_clear(d, cctx);

            /* printers: context and element strings */
            GR_MUST_SUCCEED(gr_ctx_get_str(&str, ctx));
            flint_free(str);
            GR_MUST_SUCCEED(gr_ctx_get_str(&str, cctx));
            flint_free(str);
            GR_MUST_SUCCEED(gr_get_str(&str, a, ctx));
            flint_free(str);
            c = gr_heap_init(cctx);
            GR_MUST_SUCCEED(gr_randtest(c, state, cctx));
            GR_MUST_SUCCEED(gr_get_str(&str, c, cctx));
            flint_free(str);
            gr_heap_clear(c, cctx);

            arf_clear(fa); arf_clear(fb); arf_clear(ref);
            gr_heap_clear(a, ctx);
            gr_heap_clear(b, ctx);
            gr_ctx_clear(ctx);
            gr_ctx_clear(cctx);
        }
    }

    /* the a-priori bound machinery over every plan shape */
    {
        ulong n_tab[] = { 7, 12, 15, 16, 23 };
        int alg_tab[] = { GR_DFT_ALG_NAIVE, GR_DFT_ALG_MIXED,
                GR_DFT_ALG_PFA, GR_DFT_ALG_BAILEY, GR_DFT_ALG_BLUESTEIN };
        slong bi;
        gr_ctx_t rctx, cctx;

        GR_MUST_SUCCEED(gr_dft_ctx_init_nfixed(rctx, 2));
        GR_MUST_SUCCEED(gr_dft_ctx_init_nfixed_complex(cctx, 2));

        for (bi = 0; bi < 5; bi++)
        {
            gr_dft_pre_t P;
            double peak, err;

            if (gr_dft_precomp_init_karatsuba(P, n_tab[bi], alg_tab[bi],
                    0, rctx, cctx) != GR_SUCCESS)
                TEST_FUNCTION_FAIL("bound plan init, n = %wu\n", n_tab[bi]);

            gr_dft_precomp_nfixed_bound(&peak, &err, 0.001, 1.0, P);
            if (!(peak >= 0.001) || !(err >= 1.0) || peak != peak ||
                    err != err)
                TEST_FUNCTION_FAIL("bound values, n = %wu\n", n_tab[bi]);

            gr_dft_precomp_clear(P);
        }

        gr_ctx_clear(rctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
