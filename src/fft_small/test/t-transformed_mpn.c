/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "mpn_extras.h"
#include "gr.h"
#include "fft_small.h"

/* (z, zn, sign) as produced by the conversions, as an fmpz */
static void
_get_fmpz(fmpz_t res, nn_srcptr z, slong zn, int sign)
{
    if (zn <= 0)
        fmpz_zero(res);
    else
    {
        fmpz_set_ui_array(res, z, zn);
        if (sign)
            fmpz_neg(res, res);
    }
}

static void
_rand_operand(fmpz_t f, nn_ptr a, slong * an, int * as, slong maxn,
              int allow_sign, flint_rand_t state)
{
    slong j;

    *an = 1 + n_randint(state, maxn);

    switch (n_randint(state, 8))
    {
        case 0:     /* one limb */
            *an = 1;
            a[0] = 1 + n_randint(state, 4);
            break;
        case 1:     /* zero */
            *an = 0;
            break;
        case 2:     /* all bits set */
            for (j = 0; j < *an; j++)
                a[j] = ~UWORD(0);
            break;
        case 3:     /* sparse: a single high bit */
            for (j = 0; j < *an; j++)
                a[j] = 0;
            a[*an - 1] = UWORD(1) << n_randint(state, FLINT_BITS);
            break;
        default:
            for (j = 0; j < *an; j++)
                a[j] = n_randlimb(state);
            break;
    }

    while (*an > 0 && a[*an - 1] == 0)
        (*an)--;

    *as = (allow_sign && *an > 0) ? (int) n_randint(state, 2) : 0;

    if (*an == 0)
        fmpz_zero(f);
    else
    {
        fmpz_set_ui_array(f, a, *an);
        if (*as)
            fmpz_neg(f, f);
    }
}

/*
    Exercises the transformed-integer ring against fmpz: random bilinear
    expressions of depth <= 2 with mixed signs and unbalanced operands,
    the exact and truncated conversions, and the refusal contract
    (GR_UNABLE on sign or capacity violations, GR_DOMAIN on undersized
    output buffers). Geometries are drawn to cross np = 4..8, several
    chunk sizes, and the threading thresholds of the transforms.
*/
TEST_FUNCTION_START(gr_transformed_mpn, state)
{
    for (slong iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        int large = (n_randint(state, 25) == 0);

        /* the signed reconstruction splits into per-segment windows with
           a serial stitch of the signed spills; vary the worker count so
           the boundaries land in different places */
        flint_set_num_threads(1 + n_randint(state, 4));
        int is_signed = n_randint(state, 2);
        slong terms_bound = 1 + n_randint(state, 5);
        slong opbits = large ? 20000 + n_randint(state, 200000)
                             : 64 + n_randint(state, 4000);
        slong bits_bound = 2 * opbits + 64 + FLINT_BIT_COUNT(terms_bound);
        slong maxn = opbits / FLINT_BITS;
        gr_ctx_t ctx;
        gr_ptr acc, x, y;
        fmpz_t fx, fy, ref, got;
        nn_ptr a, b, z;
        slong an, bn, zn_out, need, k, nterms;
        int as, bs, sign, status;

        if (gr_ctx_init_transformed_mpn(ctx, bits_bound, terms_bound,
                                        is_signed, 16) != GR_SUCCESS)
            continue;

        acc = gr_heap_init(ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        fmpz_init(fx);
        fmpz_init(fy);
        fmpz_init(ref);
        fmpz_init(got);
        a = flint_malloc((maxn + 1) * sizeof(ulong));
        b = flint_malloc((maxn + 1) * sizeof(ulong));

        nterms = 1 + n_randint(state, terms_bound);
        fmpz_zero(ref);
        status = GR_SUCCESS;

        for (k = 0; k < nterms && status == GR_SUCCESS; k++)
        {
            int sub = is_signed && n_randint(state, 2);

            _rand_operand(fx, a, &an, &as, maxn, is_signed, state);
            _rand_operand(fy, b, &bn, &bs, maxn, is_signed, state);

            status |= gr_transformed_mpn_set(x, a, an, as, ctx);
            status |= gr_transformed_mpn_set(y, b, bn, bs, ctx);
            if (status != GR_SUCCESS)
                break;

            if (k == 0)
            {
                status |= gr_mul(acc, x, y, ctx);
                fmpz_mul(ref, fx, fy);
            }
            else if (sub)
            {
                status |= gr_submul(acc, x, y, ctx);
                fmpz_submul(ref, fx, fy);
            }
            else
            {
                status |= gr_addmul(acc, x, y, ctx);
                fmpz_addmul(ref, fx, fy);
            }
        }

        /* an unsigned context must refuse everything that could produce a
           negative value rather than return a wrong one */
        if (status != GR_SUCCESS)
        {
            if (status != GR_UNABLE && status != GR_DOMAIN)
                TEST_FUNCTION_FAIL(
                    "unexpected status %d\n"
                    "is_signed = %d, terms_bound = %wd, bits_bound = %wd\n",
                    status, is_signed, terms_bound, bits_bound);
            goto cleanup;
        }

        need = gr_transformed_mpn_get_limbs(ctx, acc);
        z = flint_malloc(FLINT_MAX(need, 1) * sizeof(ulong));

        /* undersized output buffers must be refused, not truncated */
        if (need > 1)
        {
            slong small = n_randint(state, need);
            if (gr_transformed_mpn_get(z, small, &zn_out, &sign, acc, ctx)
                    != GR_DOMAIN && fmpz_size(ref) > small)
                TEST_FUNCTION_FAIL("undersized buffer not refused\n"
                    "need = %wd, given = %wd\n", need, small);
        }

        if (gr_transformed_mpn_get(z, need, &zn_out, &sign, acc, ctx)
                != GR_SUCCESS)
            TEST_FUNCTION_FAIL("exact conversion failed\n"
                "is_signed = %d, terms_bound = %wd, bits_bound = %wd\n",
                is_signed, terms_bound, bits_bound);

        _get_fmpz(got, z, zn_out, sign);

        if (!fmpz_equal(got, ref))
            TEST_FUNCTION_FAIL(
                "exact conversion wrong\n"
                "is_signed = %d, terms_bound = %wd, bits_bound = %wd\n"
                "ref bits = %wd, got bits = %wd\n"
                "ref sign = %d, got sign = %d\n",
                is_signed, terms_bound, bits_bound,
                (slong) fmpz_bits(ref), (slong) fmpz_bits(got),
                fmpz_sgn(ref), sign ? -1 : 1);

        /* zero must be canonical */
        if (fmpz_is_zero(ref) && (zn_out != 0 || sign != 0))
            TEST_FUNCTION_FAIL("zero not canonical: zn_out = %wd, sign = %d\n",
                zn_out, sign);

        /* exact cancellation: acc - acc = 0 */
        if (is_signed)
        {
            gr_ptr c = gr_heap_init(ctx);
            slong zn2;
            int sg2, st2;

            st2 = gr_sub(c, acc, acc, ctx);
            if (st2 == GR_SUCCESS)
            {
                slong nd2 = gr_transformed_mpn_get_limbs(ctx, c);
                nn_ptr z2 = flint_malloc(FLINT_MAX(nd2, 1) * sizeof(ulong));

                if (gr_transformed_mpn_get(z2, nd2, &zn2, &sg2, c, ctx)
                        == GR_SUCCESS && (zn2 != 0 || sg2 != 0))
                    TEST_FUNCTION_FAIL("cancellation not zero: "
                        "zn = %wd, sign = %d\n", zn2, sg2);
                flint_free(z2);
            }
            gr_heap_clear(c, ctx);
        }

        /* truncated conversion: at most one unit of error in the lowest
           returned limb, and exact for lo = 0 */
        {
            slong lo = n_randint(state, need + 2);
            slong needt = gr_transformed_mpn_get_limbs_trunc(ctx, acc, lo);
            nn_ptr zt = flint_malloc(FLINT_MAX(needt, 1) * sizeof(ulong));
            slong znt;
            int sgt;

            if (gr_transformed_mpn_get_trunc(zt, needt, &znt, &sgt, lo,
                                             acc, ctx) == GR_SUCCESS)
            {
                fmpz_t want, diff;

                fmpz_init(want);
                fmpz_init(diff);
                fmpz_abs(want, ref);
                fmpz_fdiv_q_2exp(want, want, FLINT_BITS * lo);

                _get_fmpz(got, zt, znt, 0);
                fmpz_sub(diff, got, want);
                fmpz_abs(diff, diff);

                if (fmpz_cmp_ui(diff, 1) > 0)
                    TEST_FUNCTION_FAIL(
                        "truncated conversion off by more than one unit\n"
                        "lo = %wd, need = %wd, want bits = %wd, "
                        "got bits = %wd\n",
                        lo, need, (slong) fmpz_bits(want),
                        (slong) fmpz_bits(got));

                if (lo == 0 && !fmpz_equal(got, want))
                    TEST_FUNCTION_FAIL("lo = 0 truncation not exact\n");

                if (!fmpz_is_zero(want) && !fmpz_is_zero(got) &&
                        sgt != (fmpz_sgn(ref) < 0))
                    TEST_FUNCTION_FAIL("truncated sign wrong: "
                        "got %d, ref sgn %d\n", sgt, fmpz_sgn(ref));

                fmpz_clear(want);
                fmpz_clear(diff);
            }
            flint_free(zt);
        }

        flint_free(z);

cleanup:
        gr_heap_clear(acc, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        fmpz_clear(fx);
        fmpz_clear(fy);
        fmpz_clear(ref);
        fmpz_clear(got);
        flint_free(a);
        flint_free(b);
        gr_ctx_clear(ctx);
    }

    /* Directed geometry sweep: walk the planner across chunk sizes and
       prime counts (np = 4..8) deterministically, since the per-slot
       reconstruction kernels of the signed export are selected by np and
       a wrong (np, N, M) pairing silently drops a cofactor limb at some
       geometries only. Each geometry is checked with a signed
       expression whose result is negative. */
    for (slong k = 9; k <= 21; k++)
    {
        slong opbits = WORD(1) << k;
        slong terms_bound = 4;
        slong bits_bound = 2 * opbits + 64 + FLINT_BIT_COUNT(terms_bound);
        slong maxn = opbits / FLINT_BITS;
        gr_ctx_t ctx;
        gr_ptr acc, x, y;
        fmpz_t fx, fy, ref, got;
        nn_ptr a, b, z;
        slong j, need, zn_out;
        int sign;

        if (gr_ctx_init_transformed_mpn(ctx, bits_bound, terms_bound, 1, 16)
                != GR_SUCCESS)
            continue;

        acc = gr_heap_init(ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        fmpz_init(fx);
        fmpz_init(fy);
        fmpz_init(ref);
        fmpz_init(got);
        a = flint_malloc(FLINT_MAX(maxn, 1) * sizeof(ulong));
        b = flint_malloc(FLINT_MAX(maxn, 1) * sizeof(ulong));

        for (j = 0; j < maxn; j++)
        {
            a[j] = n_randlimb(state);
            b[j] = n_randlimb(state);
        }
        a[maxn - 1] |= UWORD(1) << (FLINT_BITS - 1);
        b[maxn - 1] |= UWORD(1) << (FLINT_BITS - 1);
        fmpz_set_ui_array(fx, a, maxn);
        fmpz_set_ui_array(fy, b, maxn);

        /* x*x - y*y with |y| > |x| forces a negative result and a
           full-width cancellation through the centered lifts */
        GR_MUST_SUCCEED(gr_transformed_mpn_set(x, a, maxn, 0, ctx));
        GR_MUST_SUCCEED(gr_transformed_mpn_set(y, b, maxn, 0, ctx));
        GR_MUST_SUCCEED(gr_mul(acc, x, x, ctx));
        GR_MUST_SUCCEED(gr_submul(acc, y, y, ctx));
        fmpz_mul(ref, fx, fx);
        fmpz_submul(ref, fy, fy);

        need = gr_transformed_mpn_get_limbs(ctx, acc);
        z = flint_malloc(FLINT_MAX(need, 1) * sizeof(ulong));
        GR_MUST_SUCCEED(gr_transformed_mpn_get(z, need, &zn_out, &sign,
                                               acc, ctx));
        _get_fmpz(got, z, zn_out, sign);

        if (!fmpz_equal(got, ref))
            TEST_FUNCTION_FAIL(
                "geometry sweep failed at opbits = %wd\n"
                "bits_bound = %wd, ref bits = %wd, got bits = %wd\n",
                opbits, bits_bound, (slong) fmpz_bits(ref),
                (slong) fmpz_bits(got));

        /* the same geometry through the truncated conversion */
        {
            slong lo = need / 3;
            slong needt = gr_transformed_mpn_get_limbs_trunc(ctx, acc, lo);
            nn_ptr zt = flint_malloc(FLINT_MAX(needt, 1) * sizeof(ulong));
            slong znt;
            int sgt;
            fmpz_t want, diff;

            GR_MUST_SUCCEED(gr_transformed_mpn_get_trunc(zt, needt, &znt,
                                                    &sgt, lo, acc, ctx));
            fmpz_init(want);
            fmpz_init(diff);
            fmpz_abs(want, ref);
            fmpz_fdiv_q_2exp(want, want, FLINT_BITS * lo);
            _get_fmpz(got, zt, znt, 0);
            fmpz_sub(diff, got, want);
            fmpz_abs(diff, diff);

            if (fmpz_cmp_ui(diff, 1) > 0 || sgt != (fmpz_sgn(ref) < 0))
                TEST_FUNCTION_FAIL(
                    "truncated geometry sweep failed at opbits = %wd\n"
                    "lo = %wd, sign got %d want %d\n",
                    opbits, lo, sgt, fmpz_sgn(ref) < 0);

            fmpz_clear(want);
            fmpz_clear(diff);
            flint_free(zt);
        }

        flint_free(z);
        flint_free(a);
        flint_free(b);
        gr_heap_clear(acc, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        fmpz_clear(fx);
        fmpz_clear(fy);
        fmpz_clear(ref);
        fmpz_clear(got);
        gr_ctx_clear(ctx);
    }

    /* directed cases: exact small values, sign of the difference, and
       the refusal contract of unsigned contexts */
    {
        gr_ctx_t ctx;
        gr_ptr acc, x, y;
        ulong a[1], b[1], z[8];
        slong zn_out;
        int sign;

        if (gr_ctx_init_transformed_mpn(ctx, 4096, 4, 1, 16) == GR_SUCCESS)
        {
            acc = gr_heap_init(ctx);
            x = gr_heap_init(ctx);
            y = gr_heap_init(ctx);

            /* 1 * 1 - 1 * 2 = -1 */
            a[0] = 1;
            b[0] = 2;
            GR_MUST_SUCCEED(gr_transformed_mpn_set(x, a, 1, 0, ctx));
            GR_MUST_SUCCEED(gr_transformed_mpn_set(y, b, 1, 0, ctx));
            GR_MUST_SUCCEED(gr_mul(acc, x, x, ctx));
            GR_MUST_SUCCEED(gr_submul(acc, x, y, ctx));
            GR_MUST_SUCCEED(gr_transformed_mpn_get(z,
                    gr_transformed_mpn_get_limbs(ctx, acc),
                    &zn_out, &sign, acc, ctx));

            if (zn_out != 1 || sign != 1 || z[0] != 1)
                TEST_FUNCTION_FAIL("1*1 - 1*2 != -1: "
                    "zn = %wd, sign = %d, z0 = %wu\n", zn_out, sign, z[0]);

            gr_heap_clear(acc, ctx);
            gr_heap_clear(x, ctx);
            gr_heap_clear(y, ctx);
            gr_ctx_clear(ctx);
        }

        if (gr_ctx_init_transformed_mpn(ctx, 4096, 4, 0, 16) == GR_SUCCESS)
        {
            acc = gr_heap_init(ctx);
            x = gr_heap_init(ctx);
            y = gr_heap_init(ctx);

            a[0] = 1;
            b[0] = 2;
            GR_MUST_SUCCEED(gr_transformed_mpn_set(x, a, 1, 0, ctx));
            GR_MUST_SUCCEED(gr_transformed_mpn_set(y, b, 1, 0, ctx));

            /* negative inputs and sign-producing operations must be
               refused on an unsigned context */
            if (gr_transformed_mpn_set(x, a, 1, 1, ctx) != GR_DOMAIN)
                TEST_FUNCTION_FAIL("negative input accepted by an "
                    "unsigned context\n");

            GR_MUST_SUCCEED(gr_mul(acc, x, y, ctx));
            if (gr_submul(acc, y, y, ctx) != GR_UNABLE)
                TEST_FUNCTION_FAIL("submul accepted by an unsigned "
                    "context\n");

            gr_heap_clear(acc, ctx);
            gr_heap_clear(x, ctx);
            gr_heap_clear(y, ctx);
            gr_ctx_clear(ctx);
        }
    }

    flint_set_num_threads(1);

    /* the truncated get is documented to differ from the floor-truncated
       value by at most 1 in the lowest returned limb; drive it with
       carry-adversarial operands (all-ones, near-cancellation) as well
       as random ones */
    {
        slong iter2;

        for (iter2 = 0; iter2 < 60 * flint_test_multiplier(); iter2++)
        {
            slong n = 32 + n_randint(state, 300);
            gr_ctx_t tctx;
            gr_ptr E;
            nn_ptr a, b, c, d, z;
            fmpz_t fa, fb, fc, fd, X, Q, G;
            slong lo, zn, zl, j, e;
            int sg, sb, adv = iter2 % 3;

            if (gr_ctx_init_transformed_mpn(tctx, 128 * n + 8, 2, 1, 16)
                    != GR_SUCCESS)
                continue;

            E = flint_malloc(5 * tctx->sizeof_elem);
#define E2_(i) GR_ENTRY(E, i, tctx->sizeof_elem)
            for (j = 0; j < 5; j++)
                gr_init(E2_(j), tctx);

            a = flint_malloc(n * sizeof(ulong));
            b = flint_malloc(n * sizeof(ulong));
            c = flint_malloc(n * sizeof(ulong));
            d = flint_malloc(n * sizeof(ulong));
            for (j = 0; j < n; j++)
            {
                a[j] = n_randlimb(state);
                b[j] = n_randlimb(state);
                c[j] = n_randlimb(state);
                d[j] = n_randlimb(state);
            }
            if (adv == 1)
                for (j = 0; j < n; j++)
                    a[j] = c[j] = ~UWORD(0);
            if (adv == 2)
            {
                flint_mpn_copyi(b, a, n);
                flint_mpn_copyi(d, c, n);
                d[0] ^= 1;
            }

            GR_MUST_SUCCEED(gr_transformed_mpn_set(E2_(0), a, n, 0, tctx));
            sb = (int) n_randint(state, 2);
            GR_MUST_SUCCEED(gr_transformed_mpn_set(E2_(1), b, n, sb, tctx));
            GR_MUST_SUCCEED(gr_transformed_mpn_set(E2_(2), c, n, 0, tctx));
            GR_MUST_SUCCEED(gr_transformed_mpn_set(E2_(3), d, n, 1, tctx));
            GR_MUST_SUCCEED(gr_mul(E2_(4), E2_(0), E2_(2), tctx));
            GR_MUST_SUCCEED(gr_addmul(E2_(4), E2_(1), E2_(3), tctx));

            lo = n / 2 + n_randint(state, n);
            zn = gr_transformed_mpn_get_limbs_trunc(tctx, E2_(4), lo);
            z = flint_malloc(FLINT_MAX(zn, 1) * sizeof(ulong));
            GR_MUST_SUCCEED(gr_transformed_mpn_get_trunc(z, zn, &zl, &sg,
                    lo, E2_(4), tctx));

            fmpz_init(fa); fmpz_init(fb); fmpz_init(fc); fmpz_init(fd);
            fmpz_init(X); fmpz_init(Q); fmpz_init(G);
            fmpz_set_ui_array(fa, a, n);
            fmpz_set_ui_array(fb, b, n);
            if (sb)
                fmpz_neg(fb, fb);
            fmpz_set_ui_array(fc, c, n);
            fmpz_set_ui_array(fd, d, n);
            fmpz_mul(X, fa, fc);
            fmpz_submul(X, fb, fd);
            fmpz_fdiv_q_2exp(Q, X, FLINT_BITS * lo);

            if (zl == 0)
                fmpz_zero(G);
            else
                fmpz_set_ui_array(G, z, FLINT_ABS(zl));
            if (sg)
                fmpz_neg(G, G);
            fmpz_sub(G, G, Q);
            e = fmpz_get_si(G);
            if (FLINT_ABS(e) > 1)
                TEST_FUNCTION_FAIL("get_trunc: error %wd exceeds the "
                        "documented 1 ulp (n = %wd, lo = %wd, adv = %d)\n",
                        e, n, lo, adv);

            flint_free(a); flint_free(b); flint_free(c); flint_free(d);
            flint_free(z);
            fmpz_clear(fa); fmpz_clear(fb); fmpz_clear(fc); fmpz_clear(fd);
            fmpz_clear(X); fmpz_clear(Q); fmpz_clear(G);
            for (j = 0; j < 5; j++)
                gr_clear(E2_(j), tctx);
            flint_free(E);
            gr_ctx_clear(tctx);
        }
    }

        /* randomized ring-compliance testing: operations outside the
       representation's provisioning (depth, accumulated terms, chunk
       capacity) answer GR_UNABLE, which the framework accepts */
    {
        gr_ctx_t ctx;
        slong it;
        for (it = 0; it < 6; it++)
        {
            slong bits = 128 << n_randint(state, 6);
            slong terms = 1 + n_randint(state, 6);
            int sgn = (int) n_randint(state, 2);
            if (gr_ctx_init_transformed_mpn(ctx, bits, terms, sgn, 8)
                    == GR_SUCCESS)
            {
                {
                    /* aliased accumulation through the internal
                       temporary */
                    gr_ptr u, v;
                    ulong lu = 7, lv = 9;
                    u = gr_heap_init(ctx); v = gr_heap_init(ctx);
                    if (gr_transformed_mpn_set(u, &lu, 1, 0, ctx)
                            == GR_SUCCESS &&
                        gr_transformed_mpn_set(v, &lv, 1, sgn, ctx)
                            == GR_SUCCESS)
                    {
                        (void) gr_addmul(u, u, v, ctx);
                        (void) gr_submul(u, v, u, ctx);
                    }
                    gr_heap_clear(u, ctx); gr_heap_clear(v, ctx);
                }
                gr_test_ring(ctx, 25, 0);
                gr_ctx_clear(ctx);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
