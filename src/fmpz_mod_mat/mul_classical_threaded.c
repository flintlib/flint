/*
    Copyright (C) 2010, 2012 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_mat.h"

#define FLINT_FMPZ_MUL_CLASSICAL_CACHE_SIZE 32768 /* size of L1 cache in words */

/*
with op = 0, computes D = A*B
with op = 1, computes D = C + A*B
with op = -1, computes D = C - A*B
*/

static inline void
_fmpz_mod_mat_addmul_basic_op(fmpz ** D, fmpz ** const C, fmpz ** const A,
               fmpz ** const B, slong m, slong k, slong n, int op, fmpz_t p)
{
    slong i, j;
    fmpz_t c;

    fmpz_init(c);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            _fmpz_vec_dot_ptr(c, A[i], B, j, k);

            if (op == 1)
                fmpz_add(c, C[i] + j, c);
            else if (op == -1)
                fmpz_sub(c, C[i] + j, c);

            fmpz_mod(D[i] + j, c, p);
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
    fmpz ** A;
    fmpz ** C;
    fmpz ** D;
    fmpz * tmp;
    fmpz * p;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
    int op;
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
    fmpz ** const A = arg.A;
    fmpz ** const C = arg.C;
    fmpz ** D = arg.D;
    fmpz * tmp = arg.tmp;
    fmpz * p = arg.p;
    int op = arg.op;
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
                _fmpz_vec_dot(c, A[i], tmp + j*k, k);

                if (op == 1)
                    fmpz_add(c, C[i] + j, c);
                else if (op == -1)
                    fmpz_sub(c, C[i] + j, c);

                fmpz_mod(D[i] + j, c, p);
            }
        }
    }
}

static inline void
_fmpz_mod_mat_addmul_transpose_threaded_pool_op(fmpz ** D, fmpz ** const C,
                            fmpz ** const A, fmpz ** const B, slong m,
                                       slong k, slong n, int op, fmpz_t p,
                               thread_pool_handle * threads, slong num_threads)
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
            fmpz_set(tmp + j*k + i, B[i] + j);

    nlimbs = fmpz_size(p);

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
        args[i].C       = C;
        args[i].D       = D;
        args[i].tmp     = tmp;
        args[i].p       = p;
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
                               thread_pool_handle * threads, slong num_threads)
{
    slong m, k, n;

    m = A->mat->r;
    k = A->mat->c;
    n = B->mat->c;

    _fmpz_mod_mat_addmul_transpose_threaded_pool_op(D->mat->rows,
                (op == 0) ? NULL : C->mat->rows, A->mat->rows, B->mat->rows,
                                    m, k, n, op, D->mod, threads, num_threads);
}

void
fmpz_mod_mat_mul_classical_threaded_op(fmpz_mod_mat_t D, const fmpz_mod_mat_t C,
            const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, int op)
{
    thread_pool_handle * threads;
    slong num_threads;

    if (A->mat->c == 0)
    {
        if (op == 0)
            fmpz_mod_mat_zero(D);
        else
            fmpz_mod_mat_set(D, C);

        return;
    }

    if (A->mat->r < FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF
        || A->mat->c < FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF
        || B->mat->c < FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF)
    {
        _fmpz_mod_mat_addmul_basic_op(D->mat->rows,
                                         (op == 0) ? NULL : C->mat->rows,
                                      A->mat->rows, B->mat->rows, A->mat->r,
                                             A->mat->c, B->mat->c, op, D->mod);

        return;
    }

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    _fmpz_mod_mat_mul_classical_threaded_pool_op(D, C, A, B, op, threads, num_threads);

    flint_give_back_threads(threads, num_threads);
}

void
fmpz_mod_mat_mul_classical_threaded(fmpz_mod_mat_t C, const fmpz_mod_mat_t A,
                                        const fmpz_mod_mat_t B)
{
    fmpz_mod_mat_mul_classical_threaded_op(C, NULL, A, B, 0);
}
