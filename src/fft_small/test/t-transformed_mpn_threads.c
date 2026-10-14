/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "mpn_extras.h"
#include "thread_support.h"
#include "gr.h"
#include "fft_small.h"

/*
    The transformed-integer ring and the complex products under threads:
    a context shared by concurrent workers (both allocation strategies,
    with and without slab staging), per-worker contexts nested on one
    thread, windowed contexts with cyclic plans, and the complex products
    (on the same pool and export machinery) running alongside, with
    unrelated fft_small work interleaved. Every result is checked
    exactly, or within the documented unit for truncated windows.
*/

/* |x - y| <= 1 for x, y read as residues mod B^zn */
static int
_thr_window_close(nn_srcptr x, nn_srcptr y, slong zn)
{
    nn_ptr t = flint_malloc(zn * sizeof(ulong));
    int ok;

    mpn_sub_n(t, x, y, zn);
    if (t[zn - 1] >> (FLINT_BITS - 1))
        mpn_neg(t, t, zn);
    ok = (t[0] <= 1) && flint_mpn_zero_p(t + 1, zn - 1);
    flint_free(t);
    return ok;
}

/* floor(x / B^lo) mod B^zn */
static void
_thr_ref_window(nn_ptr z, const fmpz_t x, slong lo, slong zn)
{
    fmpz_t t, m;
    fmpz_init(t);
    fmpz_init(m);
    fmpz_fdiv_q_2exp(t, x, FLINT_BITS * lo);
    fmpz_one(m);
    fmpz_mul_2exp(m, m, FLINT_BITS * zn);
    fmpz_mod(t, t, m);
    flint_mpn_zero(z, zn);
    if (!fmpz_is_zero(t))
        fmpz_get_ui_array(z, zn, t);
    fmpz_clear(t);
    fmpz_clear(m);
}

static void
_thr_rand_operand(fmpz_t f, nn_ptr a, slong * an, int * s, slong maxn,
              flint_rand_t state)
{
    slong j;

    *an = 1 + n_randint(state, maxn);
    for (j = 0; j < *an; j++)
        a[j] = n_randlimb(state);
    if (n_randint(state, 8) == 0)
        for (j = 0; j < *an; j++)
            a[j] = ~UWORD(0);
    a[*an - 1] |= 1;
    *s = n_randint(state, 2);
    fmpz_set_ui_array(f, a, *an);
    if (*s)
        fmpz_neg(f, f);
}

typedef struct
{
    gr_ctx_struct * shared;     /* NULL: per-worker contexts */
    slong maxn;
    ulong seed;
    int engine;                 /* also run the complex products */
    int failed;
}
_thr_worker_arg;

/* a b - c d in ctx, converted out four ways and checked; returns 0 on a
   wrong value (refusals of aliased elements are checked separately) */
static int
_thr_ring_round(gr_ctx_struct * ctx, slong maxn, flint_rand_t state)
{
    slong sz = ctx->sizeof_elem, an, k, need, vn, zn, lo, w;
    gr_ptr E;
    nn_ptr t, z, r;
    fmpz_t f[4], ref, u, g;
    int s, ok = 1, sg, st;

    t = flint_malloc(maxn * sizeof(ulong));
    for (k = 0; k < 4; k++)
        fmpz_init(f[k]);
    fmpz_init(ref);
    fmpz_init(u);
    fmpz_init(g);
    GR_TMP_INIT_VEC(E, 5, ctx);

    st = GR_SUCCESS;
    for (k = 0; k < 4; k++)
    {
        _thr_rand_operand(f[k], t, &an, &s, maxn, state);
        st |= gr_transformed_mpn_set(GR_ENTRY(E, k, sz), t, an, s, ctx);
    }
    st |= gr_mul(GR_ENTRY(E, 4, sz), GR_ENTRY(E, 0, sz), GR_ENTRY(E, 1, sz), ctx);
    st |= gr_submul(GR_ENTRY(E, 4, sz), GR_ENTRY(E, 2, sz), GR_ENTRY(E, 3, sz), ctx);
    if (st != GR_SUCCESS)
    {
        ok = 0;
        goto cleanup;
    }

    fmpz_mul(ref, f[0], f[1]);
    fmpz_mul(u, f[2], f[3]);
    fmpz_sub(ref, ref, u);
    vn = fmpz_size(ref);

    /* unrelated fft_small work on this thread while elements are live */
    if (n_randint(state, 3) == 0)
    {
        slong n = 1 + n_randint(state, 4 * maxn);
        nn_ptr x = flint_malloc(n * sizeof(ulong));
        nn_ptr y = flint_malloc(2 * n * sizeof(ulong));
        flint_mpn_rrandom(x, state, n);
        x[n - 1] |= 1;
        flint_mpn_mul(y, x, n, x, n);
        flint_free(x);
        flint_free(y);
    }

    need = gr_transformed_mpn_get_limbs(ctx, GR_ENTRY(E, 4, sz));
    z = flint_malloc((need + 2) * sizeof(ulong));
    r = flint_malloc((need + 2) * sizeof(ulong));

    /* the full reconstruction width, and a window just holding the
       value (staged) */
    for (w = 0; w < 2 && ok; w++)
    {
        slong width = (w == 0) ? need : FLINT_MAX(vn, 1);
        if (gr_transformed_mpn_get(z, width, &zn, &sg, GR_ENTRY(E, 4, sz), ctx) != GR_SUCCESS)
            ok = 0;
        else
        {
            if (zn == 0)
                fmpz_zero(g);
            else
                fmpz_set_ui_array(g, z, zn);
            if (sg)
                fmpz_neg(g, g);
            ok = fmpz_equal(g, ref);
        }
    }

    /* truncated, from a random limb (staged) */
    lo = n_randint(state, vn + 1);
    if (ok && lo < vn)
    {
        slong tn = vn - lo + 1;
        if (gr_transformed_mpn_get_trunc(z, tn, &zn, &sg, lo, GR_ENTRY(E, 4, sz), ctx) != GR_SUCCESS)
            ok = 0;
        else
        {
            fmpz_abs(u, ref);
            fmpz_fdiv_q_2exp(u, u, FLINT_BITS * lo);
            if (zn == 0)
                fmpz_zero(g);
            else
                fmpz_set_ui_array(g, z, zn);
            fmpz_sub(g, g, u);
            ok = (fmpz_cmp_si(g, -1) >= 0 && fmpz_cmp_si(g, 1) <= 0);
        }
    }

    /* a two's-complement window */
    if (ok)
    {
        lo = n_randint(state, vn + 1);
        zn = 1 + n_randint(state, vn + 2 - lo);
        if (gr_transformed_mpn_get_window(z, zn, lo, GR_ENTRY(E, 4, sz), ctx) != GR_SUCCESS)
            ok = 0;
        else
        {
            _thr_ref_window(r, ref, lo, zn);
            ok = (lo == 0) ? (mpn_cmp(z, r, zn) == 0) : _thr_window_close(z, r, zn);
        }
    }

    flint_free(z);
    flint_free(r);

cleanup:
    GR_TMP_CLEAR_VEC(E, 5, ctx);
    for (k = 0; k < 4; k++)
        fmpz_clear(f[k]);
    fmpz_clear(ref);
    fmpz_clear(u);
    fmpz_clear(g);
    flint_free(t);
    return ok;
}

/* one windowed complex product by the engine, checked against the
   schoolbook product */
static int
_thr_engine_round(slong maxn, flint_rand_t state)
{
    nn_ptr a[4], zr, zi, er, ei, wr, wi;
    slong n[4], k, lo, hi, w, zrl, zil;
    int s[4], ok;

    for (k = 0; k < 4; k++)
    {
        a[k] = flint_malloc(maxn * sizeof(ulong));
        n[k] = 1 + n_randint(state, maxn);
        flint_mpn_rrandom(a[k], state, n[k]);
        s[k] = n_randint(state, 2);
        if (n_randint(state, 6) == 0)
            flint_mpn_zero(a[k], n[k]);
    }

    w = FLINT_MAX(n[0], n[1]) + FLINT_MAX(n[2], n[3]) + 1;
    er = flint_malloc(w * sizeof(ulong));
    ei = flint_malloc(w * sizeof(ulong));
    flint_mpn_mul_complex_classical(er, &zrl, ei, &zil, a[0], n[0], s[0],
            a[1], n[1], s[1], a[2], n[2], s[2], a[3], n[3], s[3]);

    lo = n_randint(state, w);
    hi = lo + 1 + n_randint(state, w - lo);
    zr = flint_malloc((hi - lo) * sizeof(ulong));
    zi = flint_malloc((hi - lo) * sizeof(ulong));
    wr = flint_malloc(w * sizeof(ulong));
    wi = flint_malloc(w * sizeof(ulong));

    /* two's complement of the exact parts, then the windows */
    flint_mpn_zero(wr, w);
    flint_mpn_zero(wi, w);
    flint_mpn_copyi(wr, er, FLINT_ABS(zrl));
    flint_mpn_copyi(wi, ei, FLINT_ABS(zil));
    if (zrl < 0)
        mpn_neg(wr, wr, w);
    if (zil < 0)
        mpn_neg(wi, wi, w);

    ok = _fft_small_mpn_mulmid_complex(zr, zi, lo, hi, a[0], n[0], s[0],
            a[1], n[1], s[1], a[2], n[2], s[2], a[3], n[3], s[3]);
    if (ok)
    {
        if (lo == 0)
            ok = (mpn_cmp(zr, wr, hi) == 0) && (mpn_cmp(zi, wi, hi) == 0);
        else
            ok = _thr_window_close(zr, wr + lo, hi - lo)
              && _thr_window_close(zi, wi + lo, hi - lo);
    }
    else
        ok = 1;     /* no plan: nothing to check */

    for (k = 0; k < 4; k++)
        flint_free(a[k]);
    flint_free(er);
    flint_free(ei);
    flint_free(zr);
    flint_free(zi);
    flint_free(wr);
    flint_free(wi);
    return ok;
}

static void
_thr_worker(slong i, void * args)
{
    _thr_worker_arg * W = ((_thr_worker_arg *) args) + i;
    flint_rand_t state;
    slong it;

    flint_rand_init(state);
    flint_rand_set_seed(state, W->seed, W->seed ^ 0x9e3779b97f4a7c15);

    for (it = 0; it < 6 && !W->failed; it++)
    {
        if (W->shared != NULL)
        {
            if (!_thr_ring_round(W->shared, W->maxn, state))
                W->failed = 1;
        }
        else
        {
            gr_ctx_t c1, c2;
            slong maxn = 1 + n_randint(state, W->maxn);
            int strat = n_randint(state, 2) ? GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER
                                            : GR_TRANSFORMED_MPN_ALLOC_MALLOC;

            if (gr_ctx_init_transformed_mpn(c1, 2 * FLINT_BITS * (maxn + 1),
                    2, 1, 5, strat) == GR_SUCCESS)
            {
                if (!_thr_ring_round(c1, maxn, state))
                    W->failed = 1;

                /* a nested context on the same thread: its reservation
                   fails and it allocates instead */
                if (n_randint(state, 2) &&
                    gr_ctx_init_transformed_mpn(c2, 2 * FLINT_BITS * (maxn + 1),
                        2, 1, 5, GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER) == GR_SUCCESS)
                {
                    if (!_thr_ring_round(c2, maxn, state) || !_thr_ring_round(c1, maxn, state))
                        W->failed = 1;
                    if (W->engine && !_thr_engine_round(maxn, state))
                        W->failed = 1;
                    gr_ctx_clear(c2);
                }
                gr_ctx_clear(c1);
            }
        }

        if (W->engine && !_thr_engine_round(W->maxn, state))
            W->failed = 1;
    }

    flint_rand_clear(state);
}

TEST_FUNCTION_START(gr_transformed_mpn_threads, state)
{
    slong iter, refusals = 0, windowed = 0;

    /* windowed contexts: products whose slots alias in a cyclic plan
       convert through windows inside the context's window only */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong n = 300 + n_randint(state, 1500), lo, hi, k, sz, an, zn;
        gr_ptr E;
        nn_ptr t, z, r;
        fmpz_t f[2], ref;
        int s;

        /* n x n products, limbs [n, n + n/4] or so: the ring's operand
           capacity is symmetric, and a window starting near the middle
           of the product leaves room for a cyclic transform */
        lo = n;
        hi = n + 1 + n_randint(state, n / 4);
        if (gr_ctx_init_transformed_mpn_window(ctx, 2 * FLINT_BITS * n,
                1, 1, 3, n_randint(state, 2) ? GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER
                                             : GR_TRANSFORMED_MPN_ALLOC_MALLOC,
                lo, hi) != GR_SUCCESS)
            continue;

        windowed++;
        sz = ctx->sizeof_elem;
        t = flint_malloc(2 * n * sizeof(ulong));
        z = flint_malloc(3 * n * sizeof(ulong));
        r = flint_malloc(3 * n * sizeof(ulong));
        fmpz_init(f[0]);
        fmpz_init(f[1]);
        fmpz_init(ref);
        GR_TMP_INIT_VEC(E, 3, ctx);

        for (k = 0; k < 2; k++)
        {
            an = n;
            flint_mpn_rrandom(t, state, an);
            t[an - 1] |= 1;
            s = n_randint(state, 2);
            fmpz_set_ui_array(f[k], t, an);
            if (s)
                fmpz_neg(f[k], f[k]);
            GR_MUST_SUCCEED(gr_transformed_mpn_set(GR_ENTRY(E, k, sz), t, an, s, ctx));
        }
        GR_MUST_SUCCEED(gr_mul(GR_ENTRY(E, 2, sz), GR_ENTRY(E, 0, sz), GR_ENTRY(E, 1, sz), ctx));
        fmpz_mul(ref, f[0], f[1]);

        /* windows inside [lo, hi) */
        {
            slong wlo = lo + n_randint(state, hi - lo);
            zn = 1 + n_randint(state, hi - wlo);
            GR_MUST_SUCCEED(gr_transformed_mpn_get_window(z, zn, wlo, GR_ENTRY(E, 2, sz), ctx));
            _thr_ref_window(r, ref, wlo, zn);
            if (!_thr_window_close(z, r, zn))
                TEST_FUNCTION_FAIL("window: n = %wd, lo = %wd, zn = %wd\n", n, wlo, zn);
        }

        /* the operands themselves never alias */
        GR_MUST_SUCCEED(gr_transformed_mpn_get_window(z, n + 1, 0, GR_ENTRY(E, 0, sz), ctx));
        _thr_ref_window(r, f[0], 0, n + 1);
        if (mpn_cmp(z, r, n + 1) != 0)
            TEST_FUNCTION_FAIL("window of an operand: n = %wd\n", n);

        /* a cyclic plan refuses conversions reading aliased slots; a
           conversion that goes through must be right */
        {
            slong lz;
            int sg, st;

            st = gr_transformed_mpn_get_window(z, 2, 0, GR_ENTRY(E, 2, sz), ctx);
            if (st == GR_SUCCESS)
            {
                _thr_ref_window(r, ref, 0, 2);
                if (mpn_cmp(z, r, 2) != 0)
                    TEST_FUNCTION_FAIL("low window: n = %wd\n", n);
            }
            else if (st == GR_UNABLE)
                refusals++;
            else
                TEST_FUNCTION_FAIL("low window status: n = %wd\n", n);

            st = gr_transformed_mpn_get_trunc(z, 2 * n + 1, &lz, &sg, 0, GR_ENTRY(E, 2, sz), ctx);
            if (st == GR_SUCCESS)
            {
                fmpz_t g;
                fmpz_init(g);
                if (lz == 0)
                    fmpz_zero(g);
                else
                    fmpz_set_ui_array(g, z, lz);
                if (sg)
                    fmpz_neg(g, g);
                if (!fmpz_equal(g, ref))
                    TEST_FUNCTION_FAIL("full conversion: n = %wd\n", n);
                fmpz_clear(g);
            }
            else if (st != GR_UNABLE)
                TEST_FUNCTION_FAIL("full conversion status: n = %wd\n", n);
        }

        GR_TMP_CLEAR_VEC(E, 3, ctx);
        fmpz_clear(f[0]);
        fmpz_clear(f[1]);
        fmpz_clear(ref);
        flint_free(t);
        flint_free(z);
        flint_free(r);
        gr_ctx_clear(ctx);
    }

    /* the cyclic plans were reached */
    if (windowed >= 10 && refusals == 0)
        TEST_FUNCTION_FAIL("no windowed context used a cyclic plan\n");

    /* concurrent workers */
    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        slong nw = 2 + n_randint(state, 6), i, maxn = 1 + n_randint(state, 500);
        _thr_worker_arg * W = flint_malloc(nw * sizeof(_thr_worker_arg));
        gr_ctx_t shared;
        int mode = n_randint(state, 3), have_shared = 0;

        flint_set_num_threads(1 + n_randint(state, 4));

        if (mode < 2)
        {
            int strat = (mode == 0) ? GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER
                                    : GR_TRANSFORMED_MPN_ALLOC_MALLOC;
            if (n_randint(state, 2))
                strat |= GR_TRANSFORMED_MPN_SCRATCH_FROM_SLAB;
            have_shared = (gr_ctx_init_transformed_mpn(shared,
                    2 * FLINT_BITS * (maxn + 1), 2, 1, 5 * nw, strat) == GR_SUCCESS);
        }

        for (i = 0; i < nw; i++)
        {
            W[i].shared = have_shared ? shared : NULL;
            W[i].maxn = maxn;
            W[i].seed = n_randlimb(state);
            W[i].engine = n_randint(state, 2);
            W[i].failed = 0;
        }

        flint_parallel_do(_thr_worker, W, nw, 0, FLINT_PARALLEL_STRIDED);

        for (i = 0; i < nw; i++)
            if (W[i].failed)
                TEST_FUNCTION_FAIL("worker %wd of %wd failed (mode %d, maxn %wd)\n",
                                   i, nw, mode, maxn);

        if (have_shared)
            gr_ctx_clear(shared);
        flint_free(W);
    }

    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
