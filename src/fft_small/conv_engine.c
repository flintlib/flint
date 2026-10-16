/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "nmod.h"
#include "fft_small.h"

/*
    The two-stage convolution engine shared by the fused multiplication
    drivers. Stage 1 computes, for each of the np primes of the plan, the
    (possibly truncated) convolution of the images of the operands: convert
    the operands into per-prime double buffers, forward transform, multiply
    pointwise with the normalization factor P->m[i] folded in, and inverse
    transform in place. Stage 2 chinese remainders BLK_SZ-aligned ranges of
    output coefficients into the destination.

    Both stages parallelize as the original per-type drivers did: stage 1
    over the primes (with each worker optionally requesting one extra
    thread to transform b concurrently with a), stage 2 over output
    ranges, with an optional re-request of up to 8 threads in between.
    The thresholds live in the caller-provided tuning table so that each
    coefficient type keeps its own tuned values.
*/

#define MAX_CONV_THREADS 8

typedef struct {
    const fft_small_plan_struct* P;
    const fft_small_conv_arg_struct* A;
    ulong start_pi;
    ulong stop_pi;
    double* abuf;
    double* bbuf;
    ulong ioff;
    int want_worker;
} s1worker_struct;

/* Convert `A->b` into `bbuf` and forward-transform it in place. Used as a
   helper for `s1worker_func` in order to process `bbuf` and `abuf` in
   parallel, assuming `!squaring`. */
static void extra_func(void* varg)
{
    s1worker_struct* X = (s1worker_struct*) varg;
    const fft_small_plan_struct* P = X->P;
    const fft_small_conv_arg_struct* A = X->A;
    sd_fft_ctx_struct* Q = P->ffts + X->ioff;

    A->mod_fn(X->bbuf, A->btrunc, A->b, A->bn, A->baux, Q, A->params);
    sd_fft_trunc(Q, X->bbuf, P->depth, A->btrunc, P->ztrunc);
}

/* Compute convolutions modulo the selected precomputed FFT primes.
   Specifically, for each `start_pi <= i < stop_pi`, reduce the operands
   modulo `P->ffts[i + P->offset].p` and convolve, writing the result
   in-place to `abuf + stride*i`.
   If `squaring` is true, the image of `a` is squared and `bbuf` is not
   used. If `A->bfft` is not NULL, it holds forward transforms of the
   second operand which are used directly.
   If `flint_request_threads` returns the requested number of threads,
   then each `s1worker_func` processes exactly one prime.  */
static void s1worker_func(void* varg)
{
    s1worker_struct* X = (s1worker_struct*) varg;
    const fft_small_plan_struct* P = X->P;
    const fft_small_conv_arg_struct* A = X->A;
    ulong i;
    thread_pool_handle* handles = NULL;
    slong nworkers = 0;

    if (X->want_worker)
        nworkers = flint_request_threads(&handles, 2);

    for (i = X->start_pi; i < X->stop_pi; i++)
    {
        ulong ioff = i + P->offset;
        double* abuf = X->abuf + P->stride*i;
        double* bbuf = X->bbuf;
        const double* bfft = NULL;
        sd_fft_ctx_struct* Q = P->ffts + ioff;

        if (A->bfft != NULL)
        {
            bfft = A->bfft + A->bfft_stride*i;
        }
        else if (!A->squaring)
        {
            if (nworkers > 0)
            {
                X->ioff = ioff; /* read by extra_func on the worker thread */
                thread_pool_wake(global_thread_pool, handles[0], 0, extra_func, X);
            }
            else
            {
                A->mod_fn(bbuf, A->btrunc, A->b, A->bn, A->baux, Q, A->params);
                sd_fft_trunc(Q, bbuf, P->depth, A->btrunc, P->ztrunc);
            }

            bfft = bbuf;
        }

        A->mod_fn(abuf, A->atrunc, A->a, A->an, A->aaux, Q, A->params);
        sd_fft_trunc(Q, abuf, P->depth, A->atrunc, P->ztrunc);

        if (A->bfft == NULL && !A->squaring && nworkers > 0)
            thread_pool_wait(global_thread_pool, handles[0]);

        if (A->squaring)
            sd_fft_ctx_point_sqr(Q, abuf, P->m[i], P->depth);
        else
            sd_fft_ctx_point_mul(Q, abuf, bfft, P->m[i], P->depth);

        sd_ifft_trunc(Q, abuf, P->depth, P->ztrunc);
    }

    flint_give_back_threads(handles, nworkers);
}

typedef struct {
    const fft_small_plan_struct* P;
    const fft_small_conv_arg_struct* A;
    void* z;
    ulong start_zi;
    ulong stop_zi;
    double* buf;
    ulong* local;
} s2worker_struct;

/* Chinese remainder an output range, writing to `z[zi - zl]` for
   `start_zi <= zi < stop_zi`.  */
static void s2worker_func(void* varg)
{
    s2worker_struct* X = (s2worker_struct*) varg;
    const fft_small_plan_struct* P = X->P;
    const fft_small_conv_arg_struct* A = X->A;

    A->crt_fn(X->z, P->zl, X->start_zi, X->stop_zi, P->ffts + P->offset,
              X->buf, P->stride, P->crts + P->offset, A->min_an_bn,
              X->local, A->params);
}

void _fft_small_conv(void* z, const fft_small_plan_t P,
                    const fft_small_conv_arg_struct* A)
{
    ulong np = P->np;
    ulong stride = P->stride;
    ulong i, o, want_threads;
    double* buf;
    thread_pool_handle* handles;
    slong nworkers;
    ulong nthreads;

    FLINT_ASSERT(P->zl < P->zh);
    FLINT_ASSERT(P->zh <= P->zn);
    FLINT_ASSERT(A->atrunc <= P->ztrunc);
    FLINT_ASSERT(A->bfft != NULL || A->btrunc <= P->ztrunc);

    want_threads = A->tuning->s1_threads(np, A->tune_n);

    nworkers = flint_request_threads(&handles, want_threads);
    nthreads = nworkers + 1;
    FLINT_ASSERT(nthreads <= MAX_CONV_THREADS);

    if (A->bfft != NULL)
        buf = (double*) mpn_ctx_fit_buffer(P->R, np*stride*sizeof(double));
    else
        buf = (double*) mpn_ctx_fit_buffer(P->R,
                                    (np + nthreads)*stride*sizeof(double));

    s1worker_struct s1args[MAX_CONV_THREADS];
    for (i = 0; i < nthreads; i++)
    {
        s1worker_struct* X = s1args + i;
        X->P = P;
        X->A = A;
        X->start_pi = (i+0)*np/nthreads;
        X->stop_pi  = (i+1)*np/nthreads;
        X->abuf = buf;
        X->bbuf = buf + (np+i)*stride;
        X->want_worker = (A->bfft == NULL) &&
                    A->tuning->s1_b_worker(np, A->tune_n, A->squaring);
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s1worker_func, s1args + i);
    s1worker_func(s1args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    if (A->tuning->s2_rethread != NULL &&
        A->tuning->s2_rethread(np, P->zn))
    {
        flint_give_back_threads(handles, nworkers);
        nworkers = flint_request_threads(&handles, MAX_CONV_THREADS);
        nthreads = nworkers + 1;
    }

    s2worker_struct s2args[MAX_CONV_THREADS];
    fft_small_crt_range_struct ranges[MAX_CONV_THREADS];
    o = P->zl;
    for (i = 0; i < nthreads; i++)
    {
        s2worker_struct* X = s2args + i;
        X->P = P;
        X->A = A;
        X->z = z;
        X->start_zi = o;
        ulong newo = n_round_down(P->zl + (i+1)*(P->zh - P->zl)/nthreads, BLK_SZ);
        o = i+1 < nthreads ? FLINT_MAX(o, newo) : P->zh;
        X->stop_zi = o;
        X->buf = buf;
        ranges[i].stop_zi = o;
        for (ulong k = 0; k < FFT_SMALL_CRT_LOCAL; k++)
            ranges[i].local[k] = 0;
        X->local = ranges[i].local;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s2worker_func, s2args + i);
    s2worker_func(s2args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);

    if (A->s2_finish != NULL)
        A->s2_finish(z, P->zl, P->zh, ranges, nthreads, A->params);
}
