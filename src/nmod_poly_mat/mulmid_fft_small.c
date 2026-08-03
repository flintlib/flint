/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "thread_pool.h"

#include <string.h>
#include <stdio.h>
#include "thread_pool.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "gr.h"
#include "gr_poly.h"

#if FLINT_HAVE_FFT_SMALL

/*
    Windowed matrix product over a ring of transformed polynomials: the
    ar*k + k*bc input entries are each transformed once and reused across
    the ar*bc*k pointwise products, and only the ar*bc output entries are
    inverse-transformed, with chinese remaindering restricted to the
    window [zl, zh). The batch is homogeneous, so the ring's shared
    transform size is tight. The ring is thread safe: the conversions in
    and the output entries are distributed over the global thread pool,
    with one accumulator element per thread.
*/

typedef struct
{
    gr_ctx_struct * ctx;
    gr_ctx_struct * tctx;
    gr_ptr TE;
    const nmod_poly_mat_struct * A;
    const nmod_poly_mat_struct * B;
    nmod_poly_mat_struct * C;
    slong ar, k, bc, nA, nB;
    slong zl, hi, wlen;
    slong tid, nthreads;
    int share;
    int status;
}
_mm_arg;

#define _MM_TE(a, i) GR_ENTRY((a)->TE, i, (a)->tctx->sizeof_elem)

static void _mm_phase1(void * varg)
{
    _mm_arg * a = (_mm_arg *) varg;
    slong t, ntot = a->nA + a->nB;

    for (t = a->tid; t < ntot; t += a->nthreads)
    {
        const nmod_poly_struct * e;
        if (t < a->nA)
            e = nmod_poly_mat_entry(a->A, t / a->k, t % a->k);
        else
        {
            slong u = t - a->nA;
            e = nmod_poly_mat_entry(a->B, u / a->bc, u % a->bc);
        }
        a->status |= _gr_set_gr_poly(_MM_TE(a, t), e->coeffs, e->length,
                                     a->ctx, a->tctx);
    }
}

static void _mm_phase2(void * varg)
{
    _mm_arg * a = (_mm_arg *) varg;
    slong nin = a->nA + a->nB;
    gr_ptr acc = _MM_TE(a, nin + a->tid);
    slong t, l;

    for (t = a->tid; t < a->ar * a->bc; t += a->nthreads)
    {
        slong i = t / a->bc, j = t % a->bc;
        nmod_poly_struct * entry = nmod_poly_mat_entry(a->C, i, j);

        a->status |= gr_mul(acc, _MM_TE(a, i * a->k),
                            a->share ? _MM_TE(a, j)
                                     : _MM_TE(a, a->nA + j), a->tctx);
        for (l = 1; l < a->k; l++)
            a->status |= gr_addmul(acc, _MM_TE(a, i * a->k + l),
                                   a->share ? _MM_TE(a, l * a->k + j)
                                            : _MM_TE(a, a->nA + l * a->bc + j),
                                   a->tctx);

        nmod_poly_fit_length(entry, a->wlen);
        /* the accumulator is this thread's own and dead after the
           entry: convert it out in place, skipping the transform copy
           and the scratch allocation; the next entry's gr_mul rebuilds
           it as a full-write destination */
        a->status |= _gr_nmod_tpoly_get_gr_poly_window_destructive(
                entry->coeffs, acc, a->zl, a->hi, a->tctx);
        if (a->wlen > a->hi - a->zl)
            _nmod_vec_zero(entry->coeffs + (a->hi - a->zl),
                           a->wlen - (a->hi - a->zl));
        _nmod_poly_set_length(entry, a->wlen);
        _nmod_poly_normalise(entry);
    }
}

static int _mulmid_fft_small_direct(nmod_poly_mat_t C,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                    slong zl, slong zh)
{
    slong ar = A->r;
    slong k = A->c;
    slong bc = B->c;
    int share = (A == B);
    slong nB;
    slong i, j, l;
    ulong amax, bmax, zn, terms, hi_prod, lo_prod;
    nmod_t mod;
    gr_ctx_t ctx, tctx;
    gr_transformed_poly_workload_t wl;
    gr_ptr TE;
    thread_pool_handle * handles = NULL;
    slong nworkers, nthreads;
    _mm_arg args[9];
    int status = GR_SUCCESS, success;

    FLINT_ASSERT(k == B->r);
    FLINT_ASSERT(C->r == ar && C->c == bc);
    FLINT_ASSERT(zl >= 0);
    FLINT_ASSERT(zh >= 0);

    if (ar == 0 || bc == 0)
        return 1;

    if (zl >= zh)
    {
        nmod_poly_mat_zero(C);
        return 1;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, ar, bc, nmod_poly_mat_modulus(A));
        success = _mulmid_fft_small_direct(T, A, B, zl, zh);
        if (success)
            nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return success;
    }

    nmod_init(&mod, nmod_poly_mat_modulus(A));

    amax = 0;
    for (i = 0; i < ar; i++)
        for (l = 0; l < k; l++)
            amax = FLINT_MAX(amax, (ulong) nmod_poly_length(nmod_poly_mat_entry(A, i, l)));
    bmax = 0;
    for (l = 0; l < k; l++)
        for (j = 0; j < bc; j++)
            bmax = FLINT_MAX(bmax, (ulong) nmod_poly_length(nmod_poly_mat_entry(B, l, j)));

    /* every product lies beyond the window, or there are no products */
    if (amax == 0 || bmax == 0 || (ulong) zl >= amax + bmax - 1)
    {
        nmod_poly_mat_zero(C);
        return 1;
    }

    zn = amax + bmax - 1;

    /* at most k*min(amax, bmax) products accumulate per coefficient */
    umul_ppmm(hi_prod, lo_prod, (ulong) k, FLINT_MIN(amax, bmax));
    if (hi_prod != 0 || lo_prod > UWORD_MAX / 8)
        return 0;
    terms = lo_prod + 2;

    gr_ctx_init_nmod(ctx, mod.n);
    nB = share ? 0 : k * bc;
    wl->num_inputs = ar * k + nB;
    wl->num_muls = ar * bc * k;
    wl->num_outputs = ar * bc;
    /* the input transforms stay live throughout (a squaring keeps a
       single pool), plus one accumulator per thread */
    wl->num_live = ar * k + nB + flint_get_num_available_threads();
    wl->mem_limit = 0;
    wl->force = 0;

    if (gr_ctx_init_gr_poly_transformed_repr(tctx, ctx, (slong) zn,
                (slong) terms, wl) != GR_SUCCESS)
    {
        gr_ctx_clear(ctx);
        return 0;
    }

    nworkers = flint_request_threads(&handles, 8);
    nthreads = nworkers + 1;

    GR_TMP_INIT_VEC(TE, ar * k + nB + nthreads, tctx);

    for (i = 0; i < nthreads; i++)
    {
        args[i].ctx = ctx; args[i].tctx = tctx; args[i].TE = TE;
        args[i].A = A; args[i].B = B; args[i].C = C;
        args[i].ar = ar; args[i].k = k; args[i].bc = bc;
        args[i].nA = ar * k;
        args[i].nB = nB;
        args[i].share = share;
        args[i].zl = zl;
        args[i].hi = FLINT_MIN(zh, (slong) zn);
        args[i].wlen = zh - zl;
        args[i].tid = i; args[i].nthreads = nthreads;
        args[i].status = GR_SUCCESS;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _mm_phase1, args + i);
    _mm_phase1(args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _mm_phase2, args + i);
    _mm_phase2(args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);

    for (i = 0; i < nthreads; i++)
        status |= args[i].status;

    GR_TMP_CLEAR_VEC(TE, ar * k + nB + nthreads, tctx);
    gr_ctx_clear(tctx);
    gr_ctx_clear(ctx);

    return status == GR_SUCCESS;
}

int nmod_poly_mat_mulmid_fft_small(nmod_poly_mat_t C,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                    slong zl, slong zh)
{
    slong ar = A->r, bc = B->c, bs, I, J;

    if (_mulmid_fft_small_direct(C, A, B, zl, zh))
        return 1;

    /* the full product exceeded the transform-storage budget: retry in
       row and column blocks over the full inner dimension (the window
       [zl, zh) of a middle product is linear in the operands, so blocks
       compose exactly), halving the block size until it fits; the inner
       dimension itself is never split */
    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        int success;
        nmod_poly_mat_init(T, ar, bc, nmod_poly_mat_modulus(A));
        success = nmod_poly_mat_mulmid_fft_small(T, A, B, zl, zh);
        if (success)
            nmod_poly_mat_swap(C, T);
        nmod_poly_mat_clear(T);
        return success;
    }

    for (bs = (FLINT_MAX(ar, bc) + 1) / 2; bs >= 1; bs /= 2)
    {
        int feasible = 1;

        for (I = 0; feasible && I < ar; I += bs)
            for (J = 0; feasible && J < bc; J += bs)
            {
                nmod_poly_mat_t Aw, Bw, Cw;
                slong ml = FLINT_MIN(bs, ar - I);
                slong nl = FLINT_MIN(bs, bc - J);

                nmod_poly_mat_window_init(Aw, A, I, 0, I + ml, A->c);
                nmod_poly_mat_window_init(Bw, B, 0, J, B->r, J + nl);
                nmod_poly_mat_window_init(Cw, C, I, J, I + ml, J + nl);
                feasible = _mulmid_fft_small_direct(Cw, Aw, Bw, zl, zh);
                nmod_poly_mat_window_clear(Aw);
                nmod_poly_mat_window_clear(Bw);
                nmod_poly_mat_window_clear(Cw);
            }

        if (feasible)
            return 1;
        if (bs == 1)
            break;
    }
    return 0;
}

#else

int nmod_poly_mat_mulmid_fft_small(nmod_poly_mat_t C,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                    slong zl, slong zh)
{
    return 0;
}

#endif
