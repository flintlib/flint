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
#include "fft_small.h"
#include "machine_vectors.h"

void fft_small_op_init(fft_small_op_t X, const fft_small_plan_t P)
{
    X->stride = P->stride;
    X->np = P->np;
    X->offset = P->offset;
    X->depth = P->depth;
    X->itrunc = 0;
    X->domain = FFT_SMALL_OP_PRIMAL;
    X->owns_data = 1;
    X->data = flint_aligned_alloc(4096,
                 n_round_up(P->np*P->stride*sizeof(double), 4096));
}

void fft_small_op_clear(fft_small_op_t X)
{
    if (X->owns_data)
        flint_aligned_free(X->data);
}

FLINT_FORCE_INLINE int _op_compatible(const fft_small_op_struct* X,
                                      const fft_small_plan_struct* P)
{
    return X->np == P->np && X->offset == P->offset &&
           X->depth == P->depth && X->stride == P->stride;
}

/*
    Pointwise kernels. Forward transforms leave entries in (-4n, 4n) and
    sd_ifft_trunc accepts entries in that range; sd_fft_ctx_point_mul
    keeps its two chained mulmods within the (-4n^2, 4n^2) input bound by
    first multiplying by m reduced to (-n/2, n/2], producing entries in
    (-3n/2, 3n/2) (see the range documentation of mulmod in
    machine_vectors.h). The kernels below do the same, and after additive
    steps reduce back into [-n, n] with reduce_to_pm1n, so that the result
    of any sequence of pointwise operations stays a valid input both to
    further pointwise operations and to sd_ifft_trunc.
*/

/* z = a*b*m */
static void _point_mul(const sd_fft_ctx_struct* Q,
    double* z, const double* a, const double* b, ulong m_, ulong depth)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong)m_, Q->p));
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    FLINT_ASSERT(depth >= LG_BLK_SZ);
    for (ulong I = 0; I < n_pow2(depth - LG_BLK_SZ); I++)
    {
        double* zx = z + sd_fft_ctx_blk_offset(I);
        const double* ax = a + sd_fft_ctx_blk_offset(I);
        const double* bx = b + sd_fft_ctx_blk_offset(I);
        ulong j = 0; do {
            vec8d x0, x1, b0, b1;
            x0 = vec8d_load(ax+j+0);
            x1 = vec8d_load(ax+j+8);
            b0 = vec8d_load(bx+j+0);
            b1 = vec8d_load(bx+j+8);
            x0 = vec8d_mulmod(x0, m, n, ninv);
            x1 = vec8d_mulmod(x1, m, n, ninv);
            x0 = vec8d_mulmod(x0, b0, n, ninv);
            x1 = vec8d_mulmod(x1, b1, n, ninv);
            vec8d_store(zx+j+0, x0);
            vec8d_store(zx+j+8, x1);
        } while (j += 16, j < BLK_SZ);
    }
}

/* z = a*a*m */
static void _point_sqr(const sd_fft_ctx_struct* Q,
    double* z, const double* a, ulong m_, ulong depth)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong)m_, Q->p));
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    FLINT_ASSERT(depth >= LG_BLK_SZ);
    for (ulong I = 0; I < n_pow2(depth - LG_BLK_SZ); I++)
    {
        double* zx = z + sd_fft_ctx_blk_offset(I);
        const double* ax = a + sd_fft_ctx_blk_offset(I);
        ulong j = 0; do {
            vec8d x0, x1;
            x0 = vec8d_load(ax+j+0);
            x1 = vec8d_load(ax+j+8);
            x0 = vec8d_mulmod(x0, x0, n, ninv);
            x1 = vec8d_mulmod(x1, x1, n, ninv);
            x0 = vec8d_mulmod(x0, m, n, ninv);
            x1 = vec8d_mulmod(x1, m, n, ninv);
            vec8d_store(zx+j+0, x0);
            vec8d_store(zx+j+8, x1);
        } while (j += 16, j < BLK_SZ);
    }
}

/* z = z +/- a*b*m */
static void _point_addmul(const sd_fft_ctx_struct* Q,
    double* z, const double* a, const double* b, ulong m_, int subtract,
    ulong depth)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong)m_, Q->p));
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    FLINT_ASSERT(depth >= LG_BLK_SZ);
    for (ulong I = 0; I < n_pow2(depth - LG_BLK_SZ); I++)
    {
        double* zx = z + sd_fft_ctx_blk_offset(I);
        const double* ax = a + sd_fft_ctx_blk_offset(I);
        const double* bx = b + sd_fft_ctx_blk_offset(I);
        ulong j = 0; do {
            vec8d x0, x1, b0, b1, z0, z1;
            x0 = vec8d_load(ax+j+0);
            x1 = vec8d_load(ax+j+8);
            b0 = vec8d_load(bx+j+0);
            b1 = vec8d_load(bx+j+8);
            z0 = vec8d_load(zx+j+0);
            z1 = vec8d_load(zx+j+8);
            x0 = vec8d_mulmod(x0, m, n, ninv);
            x1 = vec8d_mulmod(x1, m, n, ninv);
            x0 = vec8d_mulmod(x0, b0, n, ninv);
            x1 = vec8d_mulmod(x1, b1, n, ninv);
            if (subtract)
            {
                z0 = vec8d_sub(z0, x0);
                z1 = vec8d_sub(z1, x1);
            }
            else
            {
                z0 = vec8d_add(z0, x0);
                z1 = vec8d_add(z1, x1);
            }
            z0 = vec8d_reduce_to_pm1n(z0, n, ninv);
            z1 = vec8d_reduce_to_pm1n(z1, n, ninv);
            vec8d_store(zx+j+0, z0);
            vec8d_store(zx+j+8, z1);
        } while (j += 16, j < BLK_SZ);
    }
}

/* z = a +/- b, or z = -a for b == NULL */
static void _point_add(const sd_fft_ctx_struct* Q,
    double* z, const double* a, const double* b, int subtract, ulong depth)
{
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    FLINT_ASSERT(depth >= LG_BLK_SZ);
    for (ulong I = 0; I < n_pow2(depth - LG_BLK_SZ); I++)
    {
        double* zx = z + sd_fft_ctx_blk_offset(I);
        const double* ax = a + sd_fft_ctx_blk_offset(I);
        const double* bx = (b == NULL) ? NULL : b + sd_fft_ctx_blk_offset(I);
        ulong j = 0; do {
            vec8d x0, x1;
            x0 = vec8d_load(ax+j+0);
            x1 = vec8d_load(ax+j+8);
            if (bx == NULL)
            {
                x0 = vec8d_neg(x0);
                x1 = vec8d_neg(x1);
            }
            else if (subtract)
            {
                x0 = vec8d_sub(x0, vec8d_load(bx+j+0));
                x1 = vec8d_sub(x1, vec8d_load(bx+j+8));
                x0 = vec8d_reduce_to_pm1n(x0, n, ninv);
                x1 = vec8d_reduce_to_pm1n(x1, n, ninv);
            }
            else
            {
                x0 = vec8d_add(x0, vec8d_load(bx+j+0));
                x1 = vec8d_add(x1, vec8d_load(bx+j+8));
                x0 = vec8d_reduce_to_pm1n(x0, n, ninv);
                x1 = vec8d_reduce_to_pm1n(x1, n, ninv);
            }
            vec8d_store(zx+j+0, x0);
            vec8d_store(zx+j+8, x1);
        } while (j += 16, j < BLK_SZ);
    }
}

void fft_small_op_mul(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P)
{
    ulong i;

    FLINT_ASSERT(_op_compatible(A, P) && _op_compatible(B, P) &&
                 _op_compatible(Z, P));
    FLINT_ASSERT(A->domain == FFT_SMALL_OP_PRIMAL);
    FLINT_ASSERT(B->domain == FFT_SMALL_OP_PRIMAL);

    for (i = 0; i < P->np; i++)
        _point_mul(P->ffts + P->offset + i, Z->data + P->stride*i,
                A->data + P->stride*i, B->data + P->stride*i,
                P->m[i], P->depth);

    Z->domain = FFT_SMALL_OP_PRODUCT;
}

void fft_small_op_sqr(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_plan_t P)
{
    ulong i;

    FLINT_ASSERT(_op_compatible(A, P) && _op_compatible(Z, P));
    FLINT_ASSERT(A->domain == FFT_SMALL_OP_PRIMAL);

    for (i = 0; i < P->np; i++)
        _point_sqr(P->ffts + P->offset + i, Z->data + P->stride*i,
                A->data + P->stride*i, P->m[i], P->depth);

    Z->domain = FFT_SMALL_OP_PRODUCT;
}

void fft_small_op_addmul(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P)
{
    ulong i;

    FLINT_ASSERT(_op_compatible(A, P) && _op_compatible(B, P) &&
                 _op_compatible(Z, P));
    FLINT_ASSERT(A->domain == FFT_SMALL_OP_PRIMAL);
    FLINT_ASSERT(B->domain == FFT_SMALL_OP_PRIMAL);
    FLINT_ASSERT(Z->domain == FFT_SMALL_OP_PRODUCT);

    for (i = 0; i < P->np; i++)
        _point_addmul(P->ffts + P->offset + i, Z->data + P->stride*i,
                A->data + P->stride*i, B->data + P->stride*i,
                P->m[i], 0, P->depth);
}

void fft_small_op_submul(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P)
{
    ulong i;

    FLINT_ASSERT(_op_compatible(A, P) && _op_compatible(B, P) &&
                 _op_compatible(Z, P));
    FLINT_ASSERT(A->domain == FFT_SMALL_OP_PRIMAL);
    FLINT_ASSERT(B->domain == FFT_SMALL_OP_PRIMAL);
    FLINT_ASSERT(Z->domain == FFT_SMALL_OP_PRODUCT);

    for (i = 0; i < P->np; i++)
        _point_addmul(P->ffts + P->offset + i, Z->data + P->stride*i,
                A->data + P->stride*i, B->data + P->stride*i,
                P->m[i], 1, P->depth);
}

void fft_small_op_add(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P)
{
    ulong i;

    FLINT_ASSERT(_op_compatible(A, P) && _op_compatible(B, P) &&
                 _op_compatible(Z, P));
    FLINT_ASSERT(A->domain == B->domain);

    for (i = 0; i < P->np; i++)
        _point_add(P->ffts + P->offset + i, Z->data + P->stride*i,
                A->data + P->stride*i, B->data + P->stride*i, 0, P->depth);

    Z->domain = A->domain;
}

void fft_small_op_sub(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P)
{
    ulong i;

    FLINT_ASSERT(_op_compatible(A, P) && _op_compatible(B, P) &&
                 _op_compatible(Z, P));
    FLINT_ASSERT(A->domain == B->domain);

    for (i = 0; i < P->np; i++)
        _point_add(P->ffts + P->offset + i, Z->data + P->stride*i,
                A->data + P->stride*i, B->data + P->stride*i, 1, P->depth);

    Z->domain = A->domain;
}

void fft_small_op_neg(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_plan_t P)
{
    ulong i;

    FLINT_ASSERT(_op_compatible(A, P) && _op_compatible(Z, P));

    for (i = 0; i < P->np; i++)
        _point_add(P->ffts + P->offset + i, Z->data + P->stride*i,
                A->data + P->stride*i, NULL, 0, P->depth);

    Z->domain = A->domain;
}

typedef struct {
    fft_small_op_struct* X;
    const fft_small_plan_struct* P;
    ulong start_pi;
    ulong stop_pi;
} _ifft_worker_struct;

static void _ifft_worker_func(void* varg)
{
    _ifft_worker_struct* W = (_ifft_worker_struct*) varg;
    const fft_small_plan_struct* P = W->P;

    for (ulong i = W->start_pi; i < W->stop_pi; i++)
        sd_ifft_trunc(P->ffts + P->offset + i, W->X->data + P->stride*i,
                P->depth, P->ztrunc);
}

void fft_small_ifft(fft_small_op_t X, const fft_small_plan_t P)
{
    ulong i;
    thread_pool_handle* handles = NULL;
    slong nworkers = 0;
    ulong nthreads;
    _ifft_worker_struct args[MPN_CTX_NCRTS];

    FLINT_ASSERT(_op_compatible(X, P));
    FLINT_ASSERT(X->domain == FFT_SMALL_OP_PRODUCT);

    /* thread over the primes; ztrunc ~ an + bn, so this corresponds to
       the fused drivers threading from operand length ~1500 up */
    if (P->np >= 2 && P->ztrunc >= 3072)
        nworkers = flint_request_threads(&handles, P->np);
    nthreads = FLINT_MIN((ulong)(nworkers + 1), P->np);

    for (i = 0; i < nthreads; i++)
    {
        args[i].X = X;
        args[i].P = P;
        args[i].start_pi = (i+0)*P->np/nthreads;
        args[i].stop_pi  = (i+1)*P->np/nthreads;
    }

    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _ifft_worker_func, args + i);
    _ifft_worker_func(args + 0);
    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);
}
