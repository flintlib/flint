/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "gr.h"
#include "gr_dft.h"

/* Split-radix decomposition, valid over any ring given w^(n/2) = -1.
   For a transform of length m with root v = w^rstep (rstep = n/m),
   splitting the input into the even entries, the entries 4j+1 and the
   entries 4j+3 with transforms U (length m/2), Z, Z' (length m/4):

       t   = v^k Z_k,   t' = v^(3k) Z'_k,
       X_k          = U_k       + (t + t')
       X_{k+m/2}    = U_k       - (t + t')
       X_{k+m/4}    = U_{k+m/4} + v^(m/4) (t - t')
       X_{k+3m/4}   = U_{k+m/4} - v^(m/4) (t - t')

   for 0 <= k < m/4, where v^(m/4) = w^(n/4) plays the role of -i.
   Per recursion step this costs 2 general root multiplications per k
   plus one multiplication by w^(n/4); the latter is free in complex
   mode, in which case the total number of multiplications by roots of
   unity is about (n/3) log2(n) versus (n/2) log2(n) for radix-2,
   and the number of real multiplications matches the minimal
   (Winograd, WFTA) counts for n = 4, 8, 16.

   Writes the natural-order transform of the strided input vec into the
   contiguous output res; res must not alias vec. */
/* rotation by +-i on complex-mode components (res may alias x) */
static int
_gr_dft_split_rot_i(gr_ptr res, gr_srcptr x, int plus_i, gr_ptr rtmp,
        gr_ctx_struct * rctx, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong rsz = rctx->sizeof_elem;
    gr_srcptr xr = x, xi = GR_ENTRY(x, 1, rsz);
    gr_ptr rr = res, ri = GR_ENTRY(res, 1, rsz);

    status |= gr_set(rtmp, xr, rctx);
    if (plus_i)
    {
        /* i x = -xi + xr i */
        status |= gr_neg(rr, xi, rctx);
        status |= gr_set(ri, rtmp, rctx);
    }
    else
    {
        /* -i x = xi - xr i */
        status |= gr_set(rr, xi, rctx);
        status |= gr_neg(ri, rtmp, rctx);
    }

    return status;
}

/* The combine pass over k in [lo, hi) (the full pass is [0, q)),
   shared by the serial and threaded transforms. When the complex
   Karatsuba table is not in use, the twiddle multiplications run in a
   tight loop with directly walked root-table indices (e1 = k rstep
   and e2 = 3 k rstep stay below n over the whole pass since
   q rstep = n/4, so a single conditional subtraction handles the
   wraparound of the inverse walk) and plain gr_mul, avoiding the
   per-multiplication dispatch of _gr_dft_mul_root; the rotation by
   w^(n/4) keeps the dispatch, being a free special rotation in
   complex mode. With the Karatsuba table active, the dispatched
   3-multiplication path is used throughout. */
static int
_gr_dft_split_pass(gr_ptr res, ulong rstep, ulong q, ulong lo, ulong hi,
        int inverse, gr_ptr t1, gr_ptr t2, gr_ptr t3, gr_ptr rtmp,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong k;

    if (P->stage_tab != NULL)
    {
        /* packed pair table for this size m = 4q: sequential entries
           (w^(k rstep), w^(3 k rstep mod n)); the inverse walks the
           pairs in reverse, using w^(-k r) = i w^((q-k) r) and
           w^(-3k r) = -i w^(3 (q-k) r), the quarter rotations being
           free in complex mode */
        gr_ctx_struct * rctx = P->real_ctx;
        gr_srcptr pairs = GR_ENTRY(P->stage_tab, P->n - 4 * q, sz);

        for (k = lo; k < hi; k++)
        {
            gr_ptr u1 = GR_ENTRY(res, k, sz);
            gr_ptr u2 = GR_ENTRY(res, q + k, sz);
            gr_ptr zk = GR_ENTRY(res, 2 * q + k, sz);
            gr_ptr zpk = GR_ENTRY(res, 3 * q + k, sz);

            if (!inverse)
            {
                status |= gr_mul(t1, zk, GR_ENTRY(pairs, 2 * k, sz), ctx);
                status |= gr_mul(t2, zpk, GR_ENTRY(pairs, 2 * k + 1, sz), ctx);
            }
            else if (k == 0)
            {
                status |= gr_set(t1, zk, ctx);
                status |= gr_set(t2, zpk, ctx);
            }
            else
            {
                status |= gr_mul(t1, zk,
                        GR_ENTRY(pairs, 2 * (q - k), sz), ctx);
                status |= _gr_dft_split_rot_i(t1, t1, 1, rtmp, rctx, ctx);
                status |= gr_mul(t2, zpk,
                        GR_ENTRY(pairs, 2 * (q - k) + 1, sz), ctx);
                status |= _gr_dft_split_rot_i(t2, t2, 0, rtmp, rctx, ctx);
            }

            status |= gr_add(t3, t1, t2, ctx);
            status |= gr_sub(t1, t1, t2, ctx);
            status |= _gr_dft_split_rot_i(t2, t1, inverse, rtmp, rctx, ctx);

            status |= gr_sub(zk, u1, t3, ctx);
            status |= gr_add(u1, u1, t3, ctx);
            status |= gr_sub(zpk, u2, t2, ctx);
            status |= gr_add(u2, u2, t2, ctx);
        }

        return status;
    }

    if (P->wtab == NULL)
    {
        ulong n = P->n;
        ulong s1 = rstep, s3 = 3 * rstep;
        ulong e1, e2;

        if (s3 >= n)
            s3 -= n;

        if (!inverse)
        {
            e1 = lo * rstep;
            e2 = 3 * lo * rstep;
            if (e2 >= n)
                e2 -= n;
            if (e2 >= n)
                e2 -= n;
        }
        else
        {
            /* multiply by w^(-e) = w^(n - e) */
            e1 = (lo == 0) ? 0 : n - lo * rstep;
            e2 = (3 * lo * rstep) % n;
            e2 = (e2 == 0) ? 0 : n - e2;
            s1 = n - s1;
            s3 = (s3 == 0) ? 0 : n - s3;
        }

        for (k = lo; k < hi; k++)
        {
            gr_ptr u1 = GR_ENTRY(res, k, sz);
            gr_ptr u2 = GR_ENTRY(res, q + k, sz);
            gr_ptr zk = GR_ENTRY(res, 2 * q + k, sz);
            gr_ptr zpk = GR_ENTRY(res, 3 * q + k, sz);

            status |= gr_mul(t1, zk, GR_ENTRY(P->roots, e1, sz), ctx);
            status |= gr_mul(t2, zpk, GR_ENTRY(P->roots, e2, sz), ctx);

            status |= gr_add(t3, t1, t2, ctx);
            status |= gr_sub(t1, t1, t2, ctx);
            status |= _gr_dft_mul_root(t2, t1, P->n >> 2, inverse, rtmp, P);

            status |= gr_sub(zk, u1, t3, ctx);
            status |= gr_add(u1, u1, t3, ctx);
            status |= gr_sub(zpk, u2, t2, ctx);
            status |= gr_add(u2, u2, t2, ctx);

            e1 += s1;
            if (e1 >= n)
                e1 -= n;
            e2 += s3;
            if (e2 >= n)
                e2 -= n;
        }
    }
    else
    {
        for (k = lo; k < hi; k++)
        {
            gr_ptr u1 = GR_ENTRY(res, k, sz);
            gr_ptr u2 = GR_ENTRY(res, q + k, sz);
            gr_ptr zk = GR_ENTRY(res, 2 * q + k, sz);
            gr_ptr zpk = GR_ENTRY(res, 3 * q + k, sz);

            status |= _gr_dft_mul_root(t1, zk, k * rstep, inverse, rtmp, P);
            status |= _gr_dft_mul_root(t2, zpk, 3 * k * rstep, inverse, rtmp, P);

            status |= gr_add(t3, t1, t2, ctx);
            status |= gr_sub(t1, t1, t2, ctx);
            status |= _gr_dft_mul_root(t2, t1, P->n >> 2, inverse, rtmp, P);

            status |= gr_sub(zk, u1, t3, ctx);
            status |= gr_add(u1, u1, t3, ctx);
            status |= gr_sub(zpk, u2, t2, ctx);
            status |= gr_add(u2, u2, t2, ctx);
        }
    }

    return status;
}

static int
_split(gr_ptr res, gr_srcptr vec, slong stride, ulong rstep, int depth,
        int inverse, gr_ptr t1, gr_ptr t2, gr_ptr t3, gr_ptr rtmp,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong m, q;

    if (depth == 0)
        return gr_set(res, vec, ctx);

    if (depth == 1)
    {
        status |= gr_add(res, vec, GR_ENTRY(vec, stride, sz), ctx);
        status |= gr_sub(GR_ENTRY(res, 1, sz), vec, GR_ENTRY(vec, stride, sz), ctx);
        return status;
    }

    m = UWORD(1) << depth;
    q = m >> 2;

    /* U, Z, Z' into res[0, m/2), res[m/2, 3m/4), res[3m/4, m) */
    status |= _split(res, vec, 2 * stride, 2 * rstep,
            depth - 1, inverse, t1, t2, t3, rtmp, P, ctx);
    status |= _split(GR_ENTRY(res, 2 * q, sz), GR_ENTRY(vec, stride, sz),
            4 * stride, 4 * rstep, depth - 2, inverse, t1, t2, t3, rtmp, P, ctx);
    status |= _split(GR_ENTRY(res, 3 * q, sz), GR_ENTRY(vec, 3 * stride, sz),
            4 * stride, 4 * rstep, depth - 2, inverse, t1, t2, t3, rtmp, P, ctx);

    status |= _gr_dft_split_pass(res, rstep, q, 0, q, inverse,
            t1, t2, t3, rtmp, P, ctx);

    return status;
}

int
_gr_dft_split(gr_ptr res, gr_srcptr vec, int inverse,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t1, t2, t3, rtmp = NULL;

    GR_TMP_INIT3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    status |= _split(res, vec, 1, 1, P->depth, inverse,
            t1, t2, t3, rtmp, P, ctx);

    GR_TMP_CLEAR3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    return status;
}

/* Threaded execution: the three recursive sub-transforms write
   disjoint regions of res and read disjoint strided cosets of vec, so
   they can run concurrently (fork-join, partitioning the available
   worker handles between the half-length sub-transform and the two
   quarter-length ones), and the combine pass consists of q = m/4
   independent iterations that are distributed across all workers.
   Nodes without an allocation of workers, and nodes below the
   granularity bound, run the serial code above. */

typedef struct
{
    gr_ptr res;
    gr_srcptr vec;
    slong stride;
    ulong rstep;
    int depth;
    int inverse;
    thread_pool_handle * handles;
    slong num_workers;
    slong blk;
    const gr_dft_pre_struct * P;
    gr_ctx_struct * ctx;
    int status;
}
_gr_dft_split_work_t;

static int
_split_mt(gr_ptr res, gr_srcptr vec, slong stride, ulong rstep, int depth,
        int inverse, gr_ptr t1, gr_ptr t2, gr_ptr t3, gr_ptr rtmp,
        thread_pool_handle * handles, slong num_workers, slong blk,
        const gr_dft_pre_t P, gr_ctx_t ctx);

static void _gr_dft_split_quarters_leaf(void * arg);

/* task: the two quarter-length sub-transforms, run concurrently with
   the half-length one; itself forks one quarter to a sub-worker when
   it has any */
static void
_gr_dft_split_quarters_worker(void * arg)
{
    _gr_dft_split_work_t * w = arg;
    gr_ctx_struct * ctx = w->ctx;
    const gr_dft_pre_struct * P = w->P;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong q = (UWORD(1) << w->depth) >> 2;
    gr_ptr t1, t2, t3, rtmp = NULL;

    GR_TMP_INIT3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    if (w->num_workers >= 1)
    {
        _gr_dft_split_work_t sub;
        slong wz = (w->num_workers - 1) / 2;

        sub.res = GR_ENTRY(w->res, 3 * q, sz);
        sub.vec = GR_ENTRY(w->vec, 3 * w->stride, sz);
        sub.stride = 4 * w->stride;
        sub.rstep = 4 * w->rstep;
        sub.depth = w->depth - 2;
        sub.inverse = w->inverse;
        sub.handles = w->handles + 1;
        sub.num_workers = wz;
        sub.blk = w->blk;
        sub.P = P;
        sub.ctx = ctx;
        sub.status = GR_SUCCESS;

        thread_pool_wake(global_thread_pool, w->handles[0], 0,
                _gr_dft_split_quarters_leaf, &sub);

        status |= _split_mt(GR_ENTRY(w->res, 2 * q, sz),
                GR_ENTRY(w->vec, w->stride, sz), 4 * w->stride,
                4 * w->rstep, w->depth - 2, w->inverse, t1, t2, t3, rtmp,
                w->handles + 1 + wz, w->num_workers - 1 - wz, w->blk,
                P, ctx);

        thread_pool_wait(global_thread_pool, w->handles[0]);
        status |= sub.status;
    }
    else
    {
        status |= _split_mt(GR_ENTRY(w->res, 2 * q, sz),
                GR_ENTRY(w->vec, w->stride, sz), 4 * w->stride,
                4 * w->rstep, w->depth - 2, w->inverse, t1, t2, t3, rtmp,
                NULL, 0, w->blk, P, ctx);
        status |= _split_mt(GR_ENTRY(w->res, 3 * q, sz),
                GR_ENTRY(w->vec, 3 * w->stride, sz), 4 * w->stride,
                4 * w->rstep, w->depth - 2, w->inverse, t1, t2, t3, rtmp,
                NULL, 0, w->blk, P, ctx);
    }

    GR_TMP_CLEAR3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    w->status = status;
}

static void
_gr_dft_split_quarters_leaf(void * arg)
{
    _gr_dft_split_work_t * w = arg;
    gr_ctx_struct * ctx = w->ctx;
    const gr_dft_pre_struct * P = w->P;
    gr_ptr t1, t2, t3, rtmp = NULL;

    GR_TMP_INIT3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    w->status = _split_mt(w->res, w->vec, w->stride, w->rstep, w->depth,
            w->inverse, t1, t2, t3, rtmp, w->handles, w->num_workers,
            w->blk, P, ctx);

    GR_TMP_CLEAR3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);
}

/* one chunk of the combine pass over k in [lo, hi) */
typedef struct
{
    gr_ptr res;
    ulong rstep;
    ulong q;
    ulong lo;
    ulong hi;
    int inverse;
    const gr_dft_pre_struct * P;
    gr_ctx_struct * ctx;
    int status;
}
_gr_dft_split_pass_work_t;

static void
_gr_dft_split_pass_worker(void * arg)
{
    _gr_dft_split_pass_work_t * w = arg;
    gr_ctx_struct * ctx = w->ctx;
    const gr_dft_pre_struct * P = w->P;
    int status = GR_SUCCESS;
    gr_ptr t1, t2, t3, rtmp = NULL;

    GR_TMP_INIT3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    status |= _gr_dft_split_pass(w->res, w->rstep, w->q, w->lo, w->hi,
            w->inverse, t1, t2, t3, rtmp, P, ctx);

    GR_TMP_CLEAR3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    w->status = status;
}

static int
_split_mt(gr_ptr res, gr_srcptr vec, slong stride, ulong rstep, int depth,
        int inverse, gr_ptr t1, gr_ptr t2, gr_ptr t3, gr_ptr rtmp,
        thread_pool_handle * handles, slong num_workers, slong blk,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong m, q;
    _gr_dft_split_work_t qt;
    slong wq, i, nchunks;

    if (num_workers < 1 || depth < 3 ||
        ((UWORD(1) << depth) >> 2) < (ulong) blk)
        return _split(res, vec, stride, rstep, depth, inverse,
                t1, t2, t3, rtmp, P, ctx);

    m = UWORD(1) << depth;
    q = m >> 2;
    (void) m;

    /* fork the two quarter-length sub-transforms (about half the
       recursive work) to a worker with about half of the execution
       slots (rounding the quarters' share up: with W workers there
       are W + 1 slots including this thread, and giving the quarters
       fewer than 2 slots would leave both quarters serial on one
       thread while the half side idles, which caps the speedup near
       2 at W = 3); run the half-length transform here */
    wq = (num_workers + 1) / 2;
    wq = FLINT_MAX(wq, 1);

    qt.res = res;
    qt.vec = vec;
    qt.stride = stride;
    qt.rstep = rstep;
    qt.depth = depth;
    qt.inverse = inverse;
    qt.handles = handles + 1;
    qt.num_workers = wq - 1;
    qt.blk = blk;
    qt.P = P;
    qt.ctx = ctx;
    qt.status = GR_SUCCESS;

    thread_pool_wake(global_thread_pool, handles[0], 0,
            _gr_dft_split_quarters_worker, &qt);

    status |= _split_mt(res, vec, 2 * stride, 2 * rstep, depth - 1,
            inverse, t1, t2, t3, rtmp, handles + wq, num_workers - wq,
            blk, P, ctx);

    thread_pool_wait(global_thread_pool, handles[0]);
    status |= qt.status;

    /* combine pass, distributed over all workers */
    nchunks = FLINT_MIN(num_workers + 1, (slong) (q / (ulong) blk));
    nchunks = FLINT_MAX(nchunks, 1);

    if (nchunks <= 1)
    {
        _gr_dft_split_pass_work_t pw;
        pw.res = res; pw.rstep = rstep; pw.q = q;
        pw.lo = 0; pw.hi = q;
        pw.inverse = inverse; pw.P = P; pw.ctx = ctx;
        pw.status = GR_SUCCESS;
        _gr_dft_split_pass_worker(&pw);
        status |= pw.status;
    }
    else
    {
        _gr_dft_split_pass_work_t * pw =
            flint_malloc(nchunks * sizeof(_gr_dft_split_pass_work_t));

        for (i = 0; i < nchunks; i++)
        {
            pw[i].res = res; pw[i].rstep = rstep; pw[i].q = q;
            pw[i].lo = (q * (ulong) i) / (ulong) nchunks;
            pw[i].hi = (q * (ulong) (i + 1)) / (ulong) nchunks;
            pw[i].inverse = inverse; pw[i].P = P; pw[i].ctx = ctx;
            pw[i].status = GR_SUCCESS;

            if (i < nchunks - 1)
                thread_pool_wake(global_thread_pool, handles[i], 0,
                        _gr_dft_split_pass_worker, &pw[i]);
        }

        _gr_dft_split_pass_worker(&pw[nchunks - 1]);

        for (i = 0; i < nchunks - 1; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        for (i = 0; i < nchunks; i++)
            status |= pw[i].status;

        flint_free(pw);
    }

    return status;
}

/* Top-level split-radix transform with threading when worker threads
   are available (falling back to the serial transform otherwise). */
int
_gr_dft_split_threaded(gr_ptr res, gr_srcptr vec, int inverse,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong blk, max_chunks;
    int threadsafe;
    thread_pool_handle * handles = P->threads;
    slong num_workers = P->num_threads;
    thread_pool_handle * pool_handles = NULL;
    slong pool_workers = 0;
    gr_ptr t1, t2, t3, rtmp = NULL;

    threadsafe = (gr_ctx_is_threadsafe(ctx) == T_TRUE) &&
            (P->real_ctx == NULL ||
             gr_ctx_is_threadsafe(P->real_ctx) == T_TRUE);

    if (!threadsafe)
    {
        handles = NULL;
        num_workers = 0;
    }

    blk = (P->serial_block > 0) ? P->serial_block : GR_DFT_SERIAL_BLOCK_DEFAULT;
    max_chunks = (slong) (P->n / (ulong) blk);

    if (threadsafe && handles == NULL && max_chunks >= 2 &&
        flint_get_num_threads() > 1)
    {
        pool_workers = flint_request_threads(&pool_handles,
                FLINT_MIN(max_chunks, (slong) (P->n / 2)));
        handles = pool_handles;
        num_workers = pool_workers;
    }

    if (handles == NULL || num_workers < 1 || max_chunks < 2)
    {
        if (pool_handles != NULL)
            flint_give_back_threads(pool_handles, pool_workers);
        return _gr_dft_split(res, vec, inverse, P, ctx);
    }

    GR_TMP_INIT3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    status |= _split_mt(res, vec, 1, 1, P->depth, inverse,
            t1, t2, t3, rtmp, handles, num_workers, blk, P, ctx);

    GR_TMP_CLEAR3(t1, t2, t3, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    if (pool_handles != NULL)
        flint_give_back_threads(pool_handles, pool_workers);

    return status;
}
