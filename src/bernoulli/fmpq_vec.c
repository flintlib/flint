/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "fmpq.h"
#include "bernoulli.h"

static void
bernoulli_vec_compute_one_thread(fmpq * res, slong a, slong b)
{
    slong i;
    bernoulli_rev_t iter;

    if (b <= a)
        return;

    /* Even indices */
    i = b - 1;
    i -= (i % 2);
    bernoulli_rev_init(iter, i);
    for ( ; i >= a; i -= 2)
        bernoulli_rev_next(fmpq_numref(res + i - a), fmpq_denref(res + i - a), iter);
    bernoulli_rev_clear(iter);

    /* Odd indices */
    for (i = b - 1 - (b % 2); i >= a; i -= 2)
    {
        if (i == 1)
            fmpq_set_si(res + i - a, -1, 2);
        else
            fmpq_zero(res + i - a);
    }
}

typedef struct
{
    fmpq * res;
    slong a;
    slong b;
    slong block_size;
    slong num_blocks;
}
work_chunk_t;

static void
worker(slong i, void * _work)
{
    work_chunk_t work = *((work_chunk_t *) _work);
    slong a, b;

    /* reverse strided scheduling */
    i = work.num_blocks - 1 - i;

    a = work.a + i * work.block_size;
    b = FLINT_MIN(a + work.block_size, work.b);

    bernoulli_vec_compute_one_thread(work.res + a - work.a, a, b);
}

void
bernoulli_fmpq_vec_no_cache(fmpq * res, ulong a, slong num)
{
    if (a > (UWORD(1) << 31) || num > 1000000000)
    {
        flint_throw(FLINT_ERROR, "bernoulli_fmpq_vec_no_cache: excessive input\n");
    }

    if (a == 0 && num <= 128)
    {
        arith_bernoulli_number_vec(res, num);
        return;
    }

    if (num < 200 || flint_get_num_threads() == 1)
    {
        bernoulli_vec_compute_one_thread(res, a, a + num);
    }
    else
    {
        slong num_blocks, block_size;
        work_chunk_t work;

        block_size = FLINT_MAX((a + num) / 32, 128);
        num_blocks = (num + block_size - 1) / block_size;

        work.res = res;
        work.a = a;
        work.b = a + num;
        work.block_size = block_size;
        work.num_blocks = num_blocks;

        flint_parallel_do(worker, &work, num_blocks, -1, FLINT_PARALLEL_STRIDED);
    }
}
