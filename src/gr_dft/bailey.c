/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "thread_pool.h"
#include "thread_support.h"
#include "gr.h"
#include "gr_dft.h"

/* Four-step (Bailey) algorithm. Write n = n1 * n2 with n1, n2 ~ sqrt(n)
   and view the input as an n2 x n1 matrix M[j2][j1] = x[j2 * n1 + j1].
   With j = j1 + n1 * j2 and k = k2 + n2 * k1,

       X_k = sum_{j1} w^(j1 k2) [ sum_{j2} x_{j1 + n1 j2} (w^n1)^(j2 k2) ]
                     (w^n2)^(j1 k1),

   giving: (1) length-n2 transforms along the columns (root w^n1),
   (2) pointwise twiddles w^(j1 k2), (3) length-n1 transforms along the
   rows (root w^n2), (4) a final transpose from C[k2 * n1 + k1] to
   X[k2 + n2 * k1], which is skipped in scrambled mode.

   Phases (1)-(3) consist of independent units of work (columns for
   (1), rows for (2) and (3)) and are optionally distributed over
   worker threads, either attached to the plan or requested from the
   global thread pool for the duration of the transform. Each phase is
   a single scatter/join round, so a threaded transform requires only
   three synchronization points. Threads are only used when
   gr_ctx_is_threadsafe certifies that the operations of the ring (and
   of the real ring, in complex Karatsuba mode) are safe for concurrent use of
   the same context; otherwise the transform runs serially. */

/* Relocate elements according to a transposition, using flat memory
   moves (gr elements are relocatable). */
static void
_transpose(gr_ptr x, ulong n1, ulong n2, int inverse, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    ulong n = n1 * n2, k1, k2;
    char * tmp = flint_malloc(n * sz);
    char * xp = (char *) x;

    if (!inverse)
    {
        /* tmp[k2 + n2 * k1] = x[k2 * n1 + k1] */
        for (k2 = 0; k2 < n2; k2++)
            for (k1 = 0; k1 < n1; k1++)
                memcpy(tmp + (k2 + n2 * k1) * sz,
                       xp + (k2 * n1 + k1) * sz, sz);
    }
    else
    {
        /* tmp[k2 * n1 + k1] = x[k2 + n2 * k1] */
        for (k2 = 0; k2 < n2; k2++)
            for (k1 = 0; k1 < n1; k1++)
                memcpy(tmp + (k2 * n1 + k1) * sz,
                       xp + (k2 + n2 * k1) * sz, sz);
    }

    memcpy(xp, tmp, n * sz);
    flint_free(tmp);
}

#define GR_DFT_BAILEY_COLS    0  /* column transforms (sub-plan P2) */
#define GR_DFT_BAILEY_TWIDDLE 1
#define GR_DFT_BAILEY_ROWS    2  /* row transforms (sub-plan P1) */

typedef struct
{
    gr_ptr x;
    const gr_dft_pre_struct * P;
    gr_ctx_struct * ctx;
    slong lo;
    slong hi;
    int phase;
    int inverse;
    int status;
}
_gr_dft_bailey_work_t;

static void
_gr_dft_bailey_worker(void * varg)
{
    _gr_dft_bailey_work_t * arg = (_gr_dft_bailey_work_t *) varg;
    int status = GR_SUCCESS;
    gr_ctx_struct * ctx = arg->ctx;
    const gr_dft_pre_struct * P = arg->P;
    slong sz = ctx->sizeof_elem;
    ulong n1 = P->n1;
    slong i;

    if (arg->phase == GR_DFT_BAILEY_COLS)
    {
        for (i = arg->lo; i < arg->hi; i++)
            status |= _gr_dft_ct(GR_ENTRY(arg->x, i, sz), n1,
                    arg->inverse, 0, P->P2, ctx);
    }
    else if (arg->phase == GR_DFT_BAILEY_ROWS)
    {
        for (i = arg->lo; i < arg->hi; i++)
            status |= _gr_dft_ct(GR_ENTRY(arg->x, i * n1, sz), 1,
                    arg->inverse, 0, P->P1, ctx);
    }
    else
    {
        gr_ptr t, rtmp = NULL;
        slong j1;

        GR_TMP_INIT(t, ctx);
        if (P->real_ctx != NULL)
            GR_TMP_INIT(rtmp, P->real_ctx);

        for (i = FLINT_MAX(arg->lo, 1); i < arg->hi; i++)
        {
            for (j1 = 1; j1 < (slong) n1; j1++)
            {
                gr_ptr a = GR_ENTRY(arg->x, i * n1 + j1, sz);

                /* j1 * k2 <= (n1 - 1)(n2 - 1) < n, no reduction needed */
                status |= _gr_dft_mul_root(t, a, j1 * i,
                        arg->inverse, rtmp, P);
                status |= gr_set(a, t, ctx);
            }
        }

        GR_TMP_CLEAR(t, ctx);
        if (P->real_ctx != NULL)
            GR_TMP_CLEAR(rtmp, P->real_ctx);
    }

    arg->status = status;
}

/* Run one phase, distributing the units [0, items) over the main
   thread and at most num_workers workers, with at most max_chunks
   chunks in total. */
static int
_gr_dft_bailey_phase(int phase, gr_ptr x, int inverse, slong items,
        slong max_chunks, thread_pool_handle * handles, slong num_workers,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong nchunks, i;
    _gr_dft_bailey_work_t * args;

    nchunks = FLINT_MIN(num_workers + 1, items);
    nchunks = FLINT_MIN(nchunks, max_chunks);
    nchunks = FLINT_MAX(nchunks, 1);

    args = flint_malloc(nchunks * sizeof(_gr_dft_bailey_work_t));

    for (i = 0; i < nchunks; i++)
    {
        args[i].x = x;
        args[i].P = P;
        args[i].ctx = ctx;
        args[i].lo = (items * i) / nchunks;
        args[i].hi = (items * (i + 1)) / nchunks;
        args[i].phase = phase;
        args[i].inverse = inverse;
        args[i].status = GR_SUCCESS;

        if (i < nchunks - 1)
            thread_pool_wake(global_thread_pool, handles[i], 0,
                    _gr_dft_bailey_worker, &args[i]);
        else
            _gr_dft_bailey_worker(&args[i]);
    }

    for (i = 0; i < nchunks - 1; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    for (i = 0; i < nchunks; i++)
        status |= args[i].status;

    flint_free(args);
    return status;
}

int
_gr_dft_bailey(gr_ptr x, int inverse, int scrambled,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong n1 = P->n1, n2 = P->n2;
    slong blk, max_chunks;
    int threadsafe;
    thread_pool_handle * handles = P->threads;
    slong num_workers = P->num_threads;
    thread_pool_handle * pool_handles = NULL;
    slong pool_workers = 0;

    threadsafe = (gr_ctx_is_threadsafe(ctx) == T_TRUE) &&
            (P->real_ctx == NULL ||
             gr_ctx_is_threadsafe(P->real_ctx) == T_TRUE);

    if (!threadsafe)
    {
        handles = NULL;
        num_workers = 0;
    }

    /* threading granularity: at most n / serial_block chunks per phase */
    blk = (P->serial_block > 0) ? P->serial_block : GR_DFT_SERIAL_BLOCK_DEFAULT;
    max_chunks = (slong) (P->n / (ulong) blk);

    if (threadsafe && handles == NULL && max_chunks >= 2 &&
        flint_get_num_threads() > 1)
    {
        pool_workers = flint_request_threads(&pool_handles,
                FLINT_MIN(max_chunks, (slong) FLINT_MAX(n1, n2)));
        handles = pool_handles;
        num_workers = pool_workers;
    }

    if (handles == NULL || num_workers < 1 || max_chunks < 2)
    {
        num_workers = 0;
        handles = NULL;
        max_chunks = 1;
    }

    if (!inverse)
    {
        status |= _gr_dft_bailey_phase(GR_DFT_BAILEY_COLS, x, 0, n1,
                max_chunks, handles, num_workers, P, ctx);
        status |= _gr_dft_bailey_phase(GR_DFT_BAILEY_TWIDDLE, x, 0, n2,
                max_chunks, handles, num_workers, P, ctx);
        status |= _gr_dft_bailey_phase(GR_DFT_BAILEY_ROWS, x, 0, n2,
                max_chunks, handles, num_workers, P, ctx);

        if (!scrambled)
            _transpose(x, n1, n2, 0, ctx);
    }
    else
    {
        if (!scrambled)
            _transpose(x, n1, n2, 1, ctx);

        status |= _gr_dft_bailey_phase(GR_DFT_BAILEY_ROWS, x, 1, n2,
                max_chunks, handles, num_workers, P, ctx);
        status |= _gr_dft_bailey_phase(GR_DFT_BAILEY_TWIDDLE, x, 1, n2,
                max_chunks, handles, num_workers, P, ctx);
        status |= _gr_dft_bailey_phase(GR_DFT_BAILEY_COLS, x, 1, n1,
                max_chunks, handles, num_workers, P, ctx);
    }

    if (pool_handles != NULL)
        flint_give_back_threads(pool_handles, pool_workers);

    return status;
}
