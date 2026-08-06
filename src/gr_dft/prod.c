/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* DFT on a product of cyclic groups Z/cyc[0] x ... x Z/cyc[num-1],
   with row-major indexing (the first component varies slowest), as in
   acb_dft_prod: one cyclic DFT along each axis, with no twiddle
   factors in between. Along axis a, the lines have stride equal to
   the suffix product of the later component lengths; strided lines
   are gathered into a contiguous buffer, transformed with the
   component plan, and scattered back, while the lines of the last
   axis (stride 1) are transformed in place. The lines of one axis are
   independent, which is what the threaded version parallelizes over.

   Component plans are shared between axes of equal length (Dirichlet
   groups routinely repeat component sizes). */

static void
_gr_dft_prod_zero(gr_dft_prod_pre_t P)
{
    P->n = 1;
    P->num = 0;
    P->cyc = NULL;
    P->num_plans = 0;
    P->plan_of = NULL;
    P->plans = NULL;
    P->flags = 0;
    P->serial_block = 0;
}

int
_gr_dft_prod_precomp_init_layout(gr_dft_prod_pre_t P, const ulong * cyc,
        slong num, int flags, int complex_mode)
{
    int status = GR_SUCCESS;
    slong a, b;

    _gr_dft_prod_zero(P);

    P->num = num;
    P->flags = flags;
    P->cyc = flint_malloc(FLINT_MAX(num, 1) * sizeof(ulong));
    P->plan_of = flint_malloc(FLINT_MAX(num, 1) * sizeof(slong));
    P->plans = flint_malloc(FLINT_MAX(num, 1) * sizeof(gr_dft_pre_struct));

    for (a = 0; a < num; a++)
    {
        P->cyc[a] = cyc[a];
        P->n *= cyc[a];

        P->plan_of[a] = -1;
        for (b = 0; b < a; b++)
        {
            if (cyc[b] == cyc[a])
            {
                P->plan_of[a] = P->plan_of[b];
                break;
            }
        }

        if (P->plan_of[a] == -1)
        {
            P->plan_of[a] = P->num_plans;
            status |= _gr_dft_precomp_init_layout(P->plans + P->num_plans,
                    cyc[a], GR_DFT_ALG_AUTO, flags, complex_mode);
            P->num_plans++;
        }
    }

    return status;
}

int
_gr_dft_prod_precomp_realize(gr_dft_prod_pre_t P, gr_ctx_struct * real_ctx,
        gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong k;

    for (k = 0; k < P->num_plans && status == GR_SUCCESS; k++)
        status |= _gr_dft_precomp_realize(P->plans + k, real_ctx, ctx);

    return status;
}

int
gr_dft_prod_precomp_init(gr_dft_prod_pre_t P, const ulong * cyc, slong num,
        int flags, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong a, b;

    _gr_dft_prod_zero(P);

    P->num = num;
    P->flags = flags;
    P->cyc = flint_malloc(FLINT_MAX(num, 1) * sizeof(ulong));
    P->plan_of = flint_malloc(FLINT_MAX(num, 1) * sizeof(slong));
    P->plans = flint_malloc(FLINT_MAX(num, 1) * sizeof(gr_dft_pre_struct));

    for (a = 0; a < num && status == GR_SUCCESS; a++)
    {
        P->cyc[a] = cyc[a];
        P->n *= cyc[a];

        P->plan_of[a] = -1;
        for (b = 0; b < a; b++)
        {
            if (cyc[b] == cyc[a])
            {
                P->plan_of[a] = P->plan_of[b];
                break;
            }
        }

        if (P->plan_of[a] == -1)
        {
            P->plan_of[a] = P->num_plans;
            status |= gr_dft_precomp_init(P->plans + P->num_plans,
                    cyc[a], GR_DFT_ALG_AUTO, flags, ctx);
            P->num_plans++;
        }
    }

    return status;
}

/* Variant for rings without canonical roots of unity: w is a root of
   unity of the given order, which every component length must divide;
   the component of length m uses w^(order/m). */
int
gr_dft_prod_precomp_init_root(gr_dft_prod_pre_t P, gr_srcptr w, ulong order,
        const ulong * cyc, slong num, int flags, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong a, b;
    gr_ptr wa;

    _gr_dft_prod_zero(P);

    for (a = 0; a < num; a++)
        if (cyc[a] == 0 || order % cyc[a] != 0)
            return GR_DOMAIN;

    P->num = num;
    P->flags = flags;
    P->cyc = flint_malloc(FLINT_MAX(num, 1) * sizeof(ulong));
    P->plan_of = flint_malloc(FLINT_MAX(num, 1) * sizeof(slong));
    P->plans = flint_malloc(FLINT_MAX(num, 1) * sizeof(gr_dft_pre_struct));

    GR_TMP_INIT(wa, ctx);

    for (a = 0; a < num && status == GR_SUCCESS; a++)
    {
        P->cyc[a] = cyc[a];
        P->n *= cyc[a];

        P->plan_of[a] = -1;
        for (b = 0; b < a; b++)
        {
            if (cyc[b] == cyc[a])
            {
                P->plan_of[a] = P->plan_of[b];
                break;
            }
        }

        if (P->plan_of[a] == -1)
        {
            P->plan_of[a] = P->num_plans;
            status |= gr_pow_ui(wa, w, order / cyc[a], ctx);
            if (status == GR_SUCCESS)
                status |= gr_dft_precomp_init_root(P->plans + P->num_plans,
                        wa, cyc[a], GR_DFT_ALG_AUTO, flags, ctx);
            P->num_plans++;
        }
    }

    GR_TMP_CLEAR(wa, ctx);
    return status;
}

void
gr_dft_prod_precomp_clear(gr_dft_prod_pre_t P)
{
    slong k;

    for (k = 0; k < P->num_plans; k++)
        gr_dft_precomp_clear(P->plans + k);

    flint_free(P->cyc);
    flint_free(P->plan_of);
    flint_free(P->plans);
    _gr_dft_prod_zero(P);
}

void
gr_dft_prod_precomp_set_serial_block(gr_dft_prod_pre_t P, slong serial_block)
{
    P->serial_block = serial_block;
}

/* Compose the fixed-point error bounds of the component transforms:
   the output of one axis, exact copies aside, is the input of the
   next. */
void
gr_dft_prod_precomp_nfixed_bound(double * peak, double * err_ulps,
        double in_mag, double in_err, const gr_dft_prod_pre_t P)
{
    double mag = in_mag, err = in_err, mag2, err2;
    slong a;

    for (a = 0; a < P->num; a++)
    {
        if (P->cyc[a] == 1)
            continue;

        gr_dft_precomp_nfixed_bound(&mag2, &err2, mag, err,
                P->plans + P->plan_of[a]);
        mag = mag2;
        err = err2;
    }

    *peak = mag;
    *err_ulps = err;
}

/* transform the contiguous lines [llo, lhi) of one axis; tmp is a
   line buffer of m elements for the strided (inner > 1) case */
static int
_gr_dft_prod_lines(gr_ptr x, int inverse, ulong m, ulong inner,
        ulong llo, ulong lhi, const gr_dft_pre_struct * plan, gr_ptr tmp,
        gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong L, t;

    if (inner == 1)
    {
        for (L = llo; L < lhi; L++)
        {
            gr_ptr line = GR_ENTRY(x, (slong) (L * m), sz);
            status |= _gr_dft_precomp_raw(line, line, inverse, plan, ctx);
        }
    }
    else
    {
        for (L = llo; L < lhi; L++)
        {
            ulong o = L / inner, i = L % inner;
            gr_ptr base = GR_ENTRY(x, (slong) (o * m * inner + i), sz);

            for (t = 0; t < m; t++)
                status |= gr_set(GR_ENTRY(tmp, t, sz),
                        GR_ENTRY(base, (slong) (t * inner), sz), ctx);

            status |= _gr_dft_precomp_raw(tmp, tmp, inverse, plan, ctx);

            for (t = 0; t < m; t++)
                status |= gr_set(GR_ENTRY(base, (slong) (t * inner), sz),
                        GR_ENTRY(tmp, t, sz), ctx);
        }
    }

    return status;
}

typedef struct
{
    gr_ptr x;
    int inverse;
    ulong m;
    ulong inner;
    ulong llo;
    ulong lhi;
    const gr_dft_pre_struct * plan;
    gr_ctx_struct * ctx;
    int status;
}
_gr_dft_prod_work_t;

static void
_gr_dft_prod_worker(void * arg)
{
    _gr_dft_prod_work_t * a = arg;
    gr_ctx_struct * ctx = a->ctx;
    gr_ptr tmp = NULL;

    if (a->inner > 1)
        tmp = gr_heap_init_vec(a->m, ctx);

    a->status = _gr_dft_prod_lines(a->x, a->inverse, a->m, a->inner,
            a->llo, a->lhi, a->plan, tmp, ctx);

    if (tmp != NULL)
        gr_heap_clear_vec(tmp, a->m, ctx);
}

int
_gr_dft_prod_precomp_raw(gr_ptr res, gr_srcptr vec, int inverse,
        const gr_dft_prod_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong n = P->n;
    slong a;
    ulong inner;
    thread_pool_handle * handles = NULL;
    slong num_workers = 0;
    slong blk, max_chunks;

    if (res != vec)
        status |= _gr_vec_set(res, vec, n, ctx);

    if (P->num == 0 || n <= 1)
        return status;

    blk = (P->serial_block > 0) ? P->serial_block
                                : GR_DFT_SERIAL_BLOCK_DEFAULT;
    max_chunks = (slong) (n / (ulong) blk);

    if (max_chunks >= 2 && flint_get_num_threads() > 1 &&
        gr_ctx_is_threadsafe(ctx) == T_TRUE)
        num_workers = flint_request_threads(&handles, max_chunks);

    /* suffix products give the line strides */
    inner = n;
    for (a = 0; a < P->num && status == GR_SUCCESS; a++)
    {
        ulong m = P->cyc[a];
        ulong nlines;
        const gr_dft_pre_struct * plan = P->plans + P->plan_of[a];

        inner /= m;

        if (m == 1)
            continue;

        nlines = n / m;

        if (num_workers < 1 || nlines < 2)
        {
            gr_ptr tmp = NULL;

            if (inner > 1)
                tmp = gr_heap_init_vec(m, ctx);
            status |= _gr_dft_prod_lines(res, inverse, m, inner,
                    0, nlines, plan, tmp, ctx);
            if (tmp != NULL)
                gr_heap_clear_vec(tmp, m, ctx);
        }
        else
        {
            slong i, nchunks = FLINT_MIN(num_workers + 1, (slong) nlines);
            _gr_dft_prod_work_t * args;

            args = flint_malloc(nchunks * sizeof(_gr_dft_prod_work_t));

            for (i = 0; i < nchunks; i++)
            {
                args[i].x = res;
                args[i].inverse = inverse;
                args[i].m = m;
                args[i].inner = inner;
                args[i].llo = (nlines * (ulong) i) / (ulong) nchunks;
                args[i].lhi = (nlines * (ulong) (i + 1)) / (ulong) nchunks;
                args[i].plan = plan;
                args[i].ctx = ctx;
                args[i].status = GR_SUCCESS;

                if (i < nchunks - 1)
                    thread_pool_wake(global_thread_pool, handles[i], 0,
                            _gr_dft_prod_worker, &args[i]);
            }

            _gr_dft_prod_worker(&args[nchunks - 1]);

            for (i = 0; i < nchunks - 1; i++)
                thread_pool_wait(global_thread_pool, handles[i]);

            for (i = 0; i < nchunks; i++)
                status |= args[i].status;

            flint_free(args);
        }
    }

    if (handles != NULL)
        flint_give_back_threads(handles, num_workers);

    return status;
}

int
gr_dft_prod_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_prod_pre_t P,
        gr_ctx_t ctx)
{
    return _gr_dft_prod_precomp_raw(res, vec, 0, P, ctx);
}

/* the raw inverse omits the 1/n normalization of each component;
   scale by the full 1/n once at the end (as in
   gr_dft_inverse_precomp, using scalar division so that exact
   division succeeds even in rings without an inverse of n) */
int
gr_dft_prod_inverse_precomp(gr_ptr res, gr_srcptr vec,
        const gr_dft_prod_pre_t P, gr_ctx_t ctx)
{
    int status;

    status = _gr_dft_prod_precomp_raw(res, vec, 1, P, ctx);

    if (P->n >= 2)
        status |= _gr_vec_div_scalar_ui(res, res, (slong) P->n, P->n, ctx);

    return status;
}

int
gr_dft_prod(gr_ptr res, gr_srcptr vec, const ulong * cyc, slong num,
        gr_ctx_t ctx)
{
    int status;
    gr_dft_prod_pre_t P;

    status = gr_dft_prod_precomp_init(P, cyc, num, 0, ctx);
    if (status == GR_SUCCESS)
        status = gr_dft_prod_precomp(res, vec, P, ctx);
    gr_dft_prod_precomp_clear(P);

    return status;
}

int
gr_dft_prod_inverse(gr_ptr res, gr_srcptr vec, const ulong * cyc, slong num,
        gr_ctx_t ctx)
{
    int status;
    gr_dft_prod_pre_t P;

    status = gr_dft_prod_precomp_init(P, cyc, num, 0, ctx);
    if (status == GR_SUCCESS)
        status = gr_dft_prod_inverse_precomp(res, vec, P, ctx);
    gr_dft_prod_precomp_clear(P);

    return status;
}
