/*
    Copyright (C) 2010, 2012 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

#define FLINT_MUL_CLASSICAL_CACHE_SIZE 32768 /* size of L1 cache in words */

/*
with op = 0, computes D = A*B
with op = 1, computes D = C + A*B
with op = -1, computes D = C - A*B
*/

static inline void
_nmod_mat_addmul_basic_op(nn_ptr D, slong Dstride, nn_srcptr C, slong Cstride, nn_srcptr A, slong Astride,
    nn_srcptr B, slong Bstride, slong m, slong k, slong n, int op, nmod_t mod, dot_params_t params)
{
    slong i, j, l;
    ulong c;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            NMOD_VEC_DOT(c, l, k, A[i * Astride + l], B[l * Bstride + j], mod, params);

            if (op == 1)
                c = nmod_add(C[i * Cstride + j], c, mod);
            else if (op == -1)
                c = nmod_sub(C[i * Cstride + j], c, mod);

            D[i * Dstride + j] = c;
        }
    }
}

typedef struct
{
    slong block;
    volatile slong * i;
    volatile slong * j;
    slong k;
    slong m;
    slong n;
    dot_params_t params;
    nn_srcptr A;
    slong Astride;
    nn_srcptr C;
    slong Cstride;
    nn_ptr D;
    slong Dstride;
    nn_ptr tmp;
    nmod_t mod;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
    int op;
} nmod_mat_transpose_arg_t;

void
_nmod_mat_addmul_transpose_worker(void * arg_ptr)
{
    nmod_mat_transpose_arg_t arg = *((nmod_mat_transpose_arg_t *) arg_ptr);
    slong i, j, iend, jend, jstart;
    slong block = arg.block;
    slong k = arg.k;
    slong m = arg.m;
    slong n = arg.n;
    dot_params_t params = arg.params;
    nn_srcptr A = arg.A;
    slong Astride = arg.Astride;
    nn_srcptr C = arg.C;
    slong Cstride = arg.Cstride;
    nn_ptr D = arg.D;
    slong Dstride = arg.Dstride;
    nn_ptr tmp = arg.tmp;
    nmod_t mod = arg.mod;
    int op = arg.op;
    ulong c;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
        i = *arg.i;
        j = *arg.j;
        if (j >= n)
        {
            i += block;
            *arg.i = i;
            j = 0;
        }
        *arg.j = j + block;
#if FLINT_USES_PTHREAD
       pthread_mutex_unlock(arg.mutex);
#endif

        if (i >= m)
            return;

        iend = FLINT_MIN(i + block, m);
        jend = FLINT_MIN(j + block, n);
        jstart = j;

        for ( ; i < iend; i++)
        {
            for (j = jstart ; j < jend; j++)
            {
                c = _nmod_vec_dot(A + i * Astride, tmp + j*k, k, mod, params);

                if (op == 1)
                    c = nmod_add(C[i * Cstride + j], c, mod);
                else if (op == -1)
                    c = nmod_sub(C[i * Cstride + j], c, mod);

                D[i * Dstride + j] = c;
            }
        }
    }
}

static inline void
_nmod_mat_addmul_transpose_threaded_pool_op(nn_ptr D, slong Dstride, nn_srcptr C, slong Cstride,
                            nn_srcptr A, slong Astride, nn_srcptr B, slong Bstride, slong m,
                          slong k, slong n, int op, nmod_t mod, dot_params_t params,
                               thread_pool_handle * threads, slong num_threads)
{
    nn_ptr tmp;
    slong i, j, block;
    slong shared_i = 0, shared_j = 0;
    nmod_mat_transpose_arg_t * args;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    tmp = flint_malloc(sizeof(ulong) * k * n);

    /* transpose B */
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
            tmp[j*k + i] = B[i * Bstride + j];

    /* compute optimal block width */
    block = FLINT_MAX(FLINT_MIN(m/(num_threads + 1), n/(num_threads + 1)), 1);

    while (2*block*k > FLINT_MUL_CLASSICAL_CACHE_SIZE && block > 1)
        block >>= 1;

    args = flint_malloc(sizeof(nmod_mat_transpose_arg_t) * (num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].block   = block;
        args[i].i       = &shared_i;
        args[i].j       = &shared_j;
        args[i].k       = k;
        args[i].m       = m;
        args[i].n       = n;
        args[i].params  = params;
        args[i].A       = A;
        args[i].Astride = Astride;
        args[i].C       = C;
        args[i].Cstride = Cstride;
        args[i].D       = D;
        args[i].Dstride = Dstride;
        args[i].tmp     = tmp;
        args[i].mod     = mod;
#if FLINT_USES_PTHREAD
	args[i].mutex   = &mutex;
#endif
	args[i].op      = op;
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, threads[i], 0,
                _nmod_mat_addmul_transpose_worker, &args[i]);
    }

    _nmod_mat_addmul_transpose_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, threads[i]);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    flint_free(args);

    flint_free(tmp);
}

typedef struct
{
    slong block;
    volatile slong * i;
    volatile slong * j;
    slong M;
    slong K;
    slong N;
    slong Kpack;
    nn_srcptr A;
    slong Astride;
    nn_srcptr C;
    slong Cstride;
    nn_ptr D;
    slong Dstride;
    nn_ptr tmp;
    nmod_t mod;
    ulong mask;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
    int pack;
    int pack_bits;
    int op;
} nmod_mat_packed_arg_t;

void
_nmod_mat_addmul_packed_worker(void * arg_ptr)
{
    nmod_mat_packed_arg_t arg = *((nmod_mat_packed_arg_t *) arg_ptr);
    slong i, j, k, iend, jend, jstart;
    slong block = arg.block;
    slong M = arg.M;
    slong K = arg.K;
    slong N = arg.N;
    slong Kpack = arg.Kpack;
    nn_srcptr A = arg.A;
    slong Astride = arg.Astride;
    nn_srcptr C = arg.C;
    slong Cstride = arg.Cstride;
    nn_ptr D = arg.D;
    slong Dstride = arg.Dstride;
    nn_ptr tmp = arg.tmp;
    nmod_t mod = arg.mod;
    ulong mask = arg.mask;
    int pack = arg.pack;
    int pack_bits = arg.pack_bits;
    int op = arg.op;
    ulong c, d;
    nn_ptr Tptr;
    nn_srcptr Aptr;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
        i = *arg.i;
        j = *arg.j;
        if (j >= Kpack)
        {
            i += block;
            *arg.i = i;
            j = 0;
        }
        *arg.j = j + block;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (i >= M)
            return;

        iend = FLINT_MIN(i + block, M);
        jend = FLINT_MIN(j + block, Kpack);
        jstart = j;

        /* multiply */
        for ( ; i < iend; i++)
        {
            for (j = jstart; j < jend; j++)
            {
                Aptr = A + i * Astride;
                Tptr = tmp + j * N;

                c = 0;

                /* unroll by 4 */
                for (k = 0; k + 4 <= N; k += 4)
                {
                    c += Aptr[k + 0] * Tptr[k + 0];
                    c += Aptr[k + 1] * Tptr[k + 1];
                    c += Aptr[k + 2] * Tptr[k + 2];
                    c += Aptr[k + 3] * Tptr[k + 3];
                }

                for ( ; k < N; k++)
                    c += Aptr[k] * Tptr[k];

                /* unpack and reduce */
                for (k = 0; k < pack && j * pack + k < K; k++)
                {
                    d = (c >> (k * pack_bits)) & mask;
                    NMOD_RED(d, d, mod);

                    if (op == 1)
                        d = nmod_add(C[i * Cstride + j * pack + k], d, mod);
                    else if (op == -1)
                        d = nmod_sub(C[i * Cstride + j * pack + k], d, mod);

                    D[i * Dstride + j * pack + k] = d;
                }
            }
        }
    }
}

/* Assumes nlimbs = 1  <->  params.method <= _DOT1 */
static void
_nmod_mat_addmul_packed_threaded_pool_op(nn_ptr D, slong Dstride,
      nn_srcptr C, slong Cstride, nn_srcptr A, slong Astride, nn_srcptr B, slong Bstride,
          slong M, slong N, slong K, int op, nmod_t mod,
                               thread_pool_handle * threads, slong num_threads)
{
    slong i, j, k;
    slong Kpack, block;
    int pack, pack_bits;
    ulong c, mask;
    nn_ptr tmp;
    slong shared_i = 0, shared_j = 0;
    nmod_mat_packed_arg_t * args;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    /* bound unreduced entry */
    c = N * (mod.n-1) * (mod.n-1);
    pack_bits = FLINT_BIT_COUNT(c);
    pack = FLINT_BITS / pack_bits;
    Kpack = (K + pack - 1) / pack;

    if (pack_bits == FLINT_BITS)
        mask = UWORD(-1);
    else
        mask = (UWORD(1) << pack_bits) - 1;

    tmp = _nmod_vec_init(Kpack * N);

    /* pack and transpose B */
    for (i = 0; i < Kpack; i++)
    {
        for (k = 0; k < N; k++)
        {
            c = B[k * Bstride + i * pack];

            for (j = 1; j < pack && i * pack + j < K; j++)
                c |= B[k * Bstride + i * pack + j] << (pack_bits * j);

            tmp[i * N + k] = c;
        }
    }

    /* compute optimal block width */
    block = FLINT_MAX(FLINT_MIN(M/(num_threads + 1), Kpack/(num_threads + 1)), 1);

    while (2*block*N > FLINT_MUL_CLASSICAL_CACHE_SIZE && block > 1)
        block >>= 1;

    args = flint_malloc(sizeof(nmod_mat_packed_arg_t) * (num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].block     = block;
        args[i].i         = &shared_i;
        args[i].j         = &shared_j;
        args[i].M         = M;
        args[i].K         = K;
        args[i].N         = N;
        args[i].Kpack     = Kpack;
        args[i].A         = A;
        args[i].Astride   = Astride;
        args[i].C         = C;
        args[i].Cstride   = Cstride;
        args[i].D         = D;
        args[i].Dstride   = Dstride;
        args[i].tmp       = tmp;
        args[i].mod       = mod;
        args[i].mask      = mask;
#if FLINT_USES_PTHREAD
	args[i].mutex     = &mutex;
#endif
	args[i].pack      = pack;
        args[i].pack_bits = pack_bits;
        args[i].op        = op;
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, threads[i],
             0, _nmod_mat_addmul_packed_worker, &args[i]);
    }

    _nmod_mat_addmul_packed_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, threads[i]);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    flint_free(args);

    _nmod_vec_clear(tmp);
}

void
_nmod_mat_mul_classical_threaded_pool_op(nmod_mat_t D, const nmod_mat_t C,
                            const nmod_mat_t A, const nmod_mat_t B, int op,
                               thread_pool_handle * threads, slong num_threads)
{
    slong m, k, n;
    nmod_t mod;

    mod = A->mod;
    m = A->r;
    k = A->c;
    n = B->c;

    dot_params_t params = _nmod_vec_dot_params(k, mod);

    if (params.method == _DOT0)
        return;
    if (params.method == _DOT1 && m > 10 && k > 10 && n > 10)
    {
        _nmod_mat_addmul_packed_threaded_pool_op(D->entries, D->stride,
            (op == 0) ? NULL : C->entries,
            (op == 0) ? 0 : C->stride,
            A->entries, A->stride, B->entries, B->stride, m, k, n, op, D->mod, threads, num_threads);
    }
    else
    {
        _nmod_mat_addmul_transpose_threaded_pool_op(D->entries, D->stride,
            (op == 0) ? NULL : C->entries,
            (op == 0) ? 0 : C->stride,
            A->entries, A->stride, B->entries, B->stride, m, k, n, op, D->mod, params, threads, num_threads);
    }
}

void
_nmod_mat_mul_classical_threaded_op(nmod_mat_t D, const nmod_mat_t C,
            const nmod_mat_t A, const nmod_mat_t B, int op)
{
    thread_pool_handle * threads;
    slong num_threads;

    if (A->c == 0)
    {
        if (op == 0)
            nmod_mat_zero(D);
        else
            nmod_mat_set(D, C);

        return;
    }

    if (A->r < NMOD_MAT_MUL_TRANSPOSE_CUTOFF
        || A->c < NMOD_MAT_MUL_TRANSPOSE_CUTOFF
        || B->c < NMOD_MAT_MUL_TRANSPOSE_CUTOFF)
    {
        dot_params_t params = _nmod_vec_dot_params(A->c, D->mod);

        _nmod_mat_addmul_basic_op(D->entries, D->stride,
            (op == 0) ? NULL : C->entries,
            (op == 0) ? 0 : C->stride,
            A->entries, A->stride, B->entries, B->stride, A->r, A->c, B->c, op, D->mod, params);

        return;
    }

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    _nmod_mat_mul_classical_threaded_pool_op(D, C, A, B, op, threads, num_threads);

    flint_give_back_threads(threads, num_threads);
}

void
nmod_mat_mul_classical_threaded(nmod_mat_t C, const nmod_mat_t A,
                                                            const nmod_mat_t B)
{
    _nmod_mat_mul_classical_threaded_op(C, NULL, A, B, 0);
}
