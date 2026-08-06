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
#include "fmpz.h"
#include "fmpz_vec.h"
/* For FMPZ_MAT_CRT_MIN_WORK_PER_THREAD */
#include "fmpz_mat.h"

/* To do: implement optimizations currently used in the fmpz_mat version */

typedef struct
{
    fmpz *res;
    nn_srcptr *residues;
    slong num_primes;
    fmpz_comb_t *comb;
    slong a;
    slong b;
    int sign;
} _multi_CRT_chunk_t;

static void _multi_CRT_worker(void *_work)
{
    _multi_CRT_chunk_t *work = (_multi_CRT_chunk_t *) _work;
    fmpz_comb_temp_t temp;
    slong i, k;

    ulong * col = flint_malloc(work->num_primes * sizeof(ulong));

    fmpz_comb_temp_init(temp, *work->comb);

    for (k = work->a; k < work->b; k++)
    {
        for (i = 0; i < work->num_primes; i++)
            col[i] = work->residues[i][k];

        fmpz_multi_CRT_ui(work->res + k, col, *work->comb, temp, work->sign);
    }

    fmpz_comb_temp_clear(temp);
    flint_free(col);
}

void _fmpz_vec_multi_CRT_ui(fmpz * res, nn_srcptr * residues, slong len,
                             nn_srcptr primes, slong num_primes, int sign)
{
    slong thread_limit;

    if (num_primes == 0)
    {
        _fmpz_vec_zero(res, len);
        return;
    }

    if (num_primes == 1)
    {
        slong i;

        if (sign)
            for (i = 0; i < len; i++)
                fmpz_set_ui_smod(res + i, residues[0][i], primes[0]);
        else
            for (i = 0; i < len; i++)
                fmpz_set_ui(res + i, residues[0][i]);
        return;
    }

    thread_limit = (len * num_primes) / FMPZ_MAT_CRT_MIN_WORK_PER_THREAD;

    if (thread_limit > 1)
        thread_limit = FLINT_MIN(thread_limit, flint_get_num_threads());

    fmpz_comb_t comb;
    fmpz_comb_init(comb, primes, num_primes);

    if (thread_limit <= 1)
    {
        _multi_CRT_chunk_t work;
        work.res        = res;
        work.residues   = residues;
        work.num_primes = num_primes;
        work.comb       = &comb;
        work.a          = 0;
        work.b          = len;
        work.sign       = sign;
        _multi_CRT_worker(&work);
    }
    else
    {
        thread_pool_handle *handles;
        slong num_workers = flint_request_threads(&handles, thread_limit);
        slong num_threads = num_workers + 1;
        slong chunk_size = (len + num_threads - 1) / num_threads;
        slong i;

        TMP_INIT;
        TMP_START;

        _multi_CRT_chunk_t *work = TMP_ALLOC(num_threads * sizeof(_multi_CRT_chunk_t));

        for (i = 0; i < num_threads; i++)
        {
            work[i].res        = res;
            work[i].residues   = residues;
            work[i].num_primes = num_primes;
            work[i].comb       = &comb;
            work[i].a          = i * chunk_size;
            work[i].b          = FLINT_MIN((i + 1) * chunk_size, len);
            work[i].sign       = sign;
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _multi_CRT_worker, &work[i]);

        _multi_CRT_worker(&work[num_workers]);

        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        flint_give_back_threads(handles, num_workers);
        TMP_END;
    }

    fmpz_comb_clear(comb);
}

