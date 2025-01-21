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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

#define FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF 10

#define FLINT_FMPZ_MUL_CLASSICAL_CACHE_SIZE 32768 /* size of L1 cache in words */

/*
with op = 0, computes D = A*B
with op = 1, computes D = C + A*B
with op = -1, computes D = C - A*B
*/

static inline void
_fmpz_mod_mat_addmul_basic_op(fmpz * D, slong Dstride, const fmpz * C, slong Cstride, const fmpz * A, slong Astride,
               const fmpz * B, slong Bstride, slong m, slong k, slong n, int op, const fmpz_mod_ctx_t ctx)
{
    slong i, j, l;
    fmpz_t c;

    fmpz_init(c);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpz_zero(c);
            for (l = 0; l < k; l++)
                fmpz_addmul(c, A + Astride * i + l, B + Bstride * l + j);

            if (op == 1)
                fmpz_add(c, C + i * Cstride + j, c);
            else if (op == -1)
                fmpz_sub(c, C + i * Cstride + j, c);

            fmpz_mod_set_fmpz(D + i * Dstride + j, c, ctx);
        }
    }

    fmpz_clear(c);
}

typedef struct
{
    slong block;
    volatile slong * i;
    volatile slong * j;
    slong k;
    slong m;
    slong n;
    fmpz * A;
    slong Astride;
    fmpz * C;
    slong Cstride;
    fmpz * D;
    slong Dstride;
    fmpz * tmp;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
    int op;
    const fmpz_mod_ctx_struct * ctx;
} fmpz_mod_mat_transpose_arg_t;

void
_fmpz_mod_mat_addmul_transpose_worker(void * arg_ptr)
{
    fmpz_mod_mat_transpose_arg_t arg = *((fmpz_mod_mat_transpose_arg_t *) arg_ptr);
    slong i, j, iend, jend, jstart;
    slong block = arg.block;
    slong k = arg.k;
    slong m = arg.m;
    slong n = arg.n;
    fmpz * A = arg.A;
    slong Astride = arg.Astride;
    fmpz * C = arg.C;
    slong Cstride = arg.Cstride;
    fmpz * D = arg.D;
    slong Dstride = arg.Dstride;
    fmpz * tmp = arg.tmp;
    int op = arg.op;
    const fmpz_mod_ctx_struct * ctx = arg.ctx;
    fmpz_t c;

    fmpz_init(c);

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
        {
            fmpz_clear(c);
            return;
        }

        iend = FLINT_MIN(i + block, m);
        jend = FLINT_MIN(j + block, n);
        jstart = j;

        for ( ; i < iend; i++)
        {
            for (j = jstart ; j < jend; j++)
            {
                _fmpz_vec_dot(c, A + i * Astride, tmp + j*k, k);

                /* todo: use dot_general */
                if (op == 1)
                    fmpz_add(c, C + i * Cstride + j, c);
                else if (op == -1)
                    fmpz_sub(c, C + i * Cstride + j, c);

                fmpz_mod_set_fmpz(D + i * Dstride + j, c, ctx);
            }
        }
    }
}

static inline void
_fmpz_mod_mat_addmul_transpose_threaded_pool_op(fmpz * D, slong Dstride, fmpz * C, slong Cstride,
                            fmpz * A, slong Astride, fmpz * B, slong Bstride, slong m,
                                       slong k, slong n, int op,
                               thread_pool_handle * threads, slong num_threads, const fmpz_mod_ctx_t ctx)
{
    fmpz * tmp;
    slong i, j, block, nlimbs;
    slong shared_i = 0, shared_j = 0;
    fmpz_mod_mat_transpose_arg_t * args;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    tmp = _fmpz_vec_init(k*n);

    /* transpose B */
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
            fmpz_set(tmp + j*k + i, B + i * Bstride + j);

    nlimbs = fmpz_size(ctx->n);

    /* compute optimal block width */
    block = FLINT_MAX(FLINT_MIN(m/(num_threads + 1), n/(num_threads + 1)), 1);

    while (2*block*k*nlimbs > FLINT_FMPZ_MUL_CLASSICAL_CACHE_SIZE && block > 1)
        block >>= 1;

    args = flint_malloc(sizeof(fmpz_mod_mat_transpose_arg_t) * (num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].block   = block;
        args[i].i       = &shared_i;
        args[i].j       = &shared_j;
        args[i].k       = k;
        args[i].m       = m;
        args[i].n       = n;
        args[i].A       = A;
        args[i].Astride = Astride;
        args[i].C       = C;
        args[i].Cstride = Cstride;
        args[i].D       = D;
        args[i].Dstride = Dstride;
        args[i].tmp     = tmp;
#if FLINT_USES_PTHREAD
        args[i].mutex   = &mutex;
#endif
        args[i].op      = op;
        args[i].ctx     = ctx;
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, threads[i], 0,
                _fmpz_mod_mat_addmul_transpose_worker, &args[i]);
    }

    _fmpz_mod_mat_addmul_transpose_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, threads[i]);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    flint_free(args);

    _fmpz_vec_clear(tmp, k*n);
}

void
_fmpz_mod_mat_mul_classical_threaded_pool_op(fmpz_mod_mat_t D, const fmpz_mod_mat_t C,
                            const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, int op,
                               thread_pool_handle * threads, slong num_threads, const fmpz_mod_ctx_t ctx)
{
    slong m, k, n;

    m = A->r;
    k = A->c;
    n = B->c;

    _fmpz_mod_mat_addmul_transpose_threaded_pool_op(D->entries, D->stride,
                (op == 0) ? NULL : C->entries,
                (op == 0) ? 0 : C->stride,
                A->entries, A->stride, B->entries, B->stride,
                                    m, k, n, op, threads, num_threads, ctx);
}

void
fmpz_mod_mat_mul_classical_threaded_op(fmpz_mod_mat_t D, const fmpz_mod_mat_t C,
            const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, int op, const fmpz_mod_ctx_t ctx)
{
    thread_pool_handle * threads;
    slong num_threads;

    if (A->c == 0)
    {
        if (op == 0)
            fmpz_mod_mat_zero(D, ctx);
        else
            fmpz_mod_mat_set(D, C, ctx);
        return;
    }

    if (A->r < FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF
        || A->c < FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF
        || B->c < FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF)
    {
        _fmpz_mod_mat_addmul_basic_op(D->entries, D->stride,
                (op == 0) ? NULL : C->entries,
                (op == 0) ? 0 : C->stride,
                A->entries, A->stride, B->entries, B->stride,
                    A->r, A->c, B->c, op, ctx);
        return;
    }

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    _fmpz_mod_mat_mul_classical_threaded_pool_op(D, C, A, B, op, threads, num_threads, ctx);

    flint_give_back_threads(threads, num_threads);
}

void
fmpz_mod_mat_mul_classical_threaded(fmpz_mod_mat_t C, const fmpz_mod_mat_t A,
                                        const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_mat_mul_classical_threaded_op(C, NULL, A, B, 0, ctx);
}
