/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

/*
    Usage: p-invert [min_bits [max_bits]]

    Ratios time(naive) / time(shoup) and time(naive) / time(generic) of
    _nmod_vec_invert_naive, _shoup and _generic (above 1: faster than the
    naive inversions one by one), for a random prime of each bit size and
    several lengths. Each timing is the minimum over a few batches of about
    1 ms of calls on consecutive vectors of a pool of 2^20 random nonzero
    entries, drawn once per bit size.

    The pool must be large: the inversions (extended Euclid) are full of
    data-dependent branches, and when the same few hundred inputs come back
    the branch predictor learns them, which can make an inversion 2-4x
    faster; with a small pool, the batch versions (one inversion per
    vector) would benefit from this much more than the naive one (len
    inversions per vector), giving ratios above len.
*/

#define POOL_ENTRIES (WORD(1) << 20)
#define BATCH_US 1000
#define ROUNDS 3

typedef void (*invert_t)(nn_ptr, nn_srcptr, ulong, nmod_t);

/* nanoseconds per call of f on consecutive vectors of the pool, *v being
   the index of the next vector (kept from one call to the next, so that
   the inputs do not come back soon) */
static double
time_invert(invert_t f, nn_ptr res, nn_srcptr pool, slong nvec, ulong len,
            nmod_t mod, slong * v)
{
    timeit_t t;
    slong reps, i, r;
    double best;

#define BATCH(nreps) \
    do { \
        timeit_start_us(t); \
        for (i = 0; i < (nreps); i++) \
        { \
            f(res, pool + (*v) * len, len, mod); \
            *v = (*v + 1 == nvec) ? 0 : *v + 1; \
        } \
        timeit_stop_us(t); \
    } while (0)

    /* number of calls for a batch of at least BATCH_US microseconds */
    for (reps = 1; ; reps *= 2)
    {
        BATCH(reps);
        if (t->wall >= BATCH_US)
            break;
    }

    best = (double) t->wall / reps;
    for (r = 1; r < ROUNDS; r++)
    {
        BATCH(reps);
        if ((double) t->wall / reps < best)
            best = (double) t->wall / reps;
    }

#undef BATCH

    return 1000.0 * best;
}

int main(int argc, char ** argv)
{
    const ulong lens[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                          16, 1024, 65536};
    const int nlens = sizeof(lens) / sizeof(lens[0]);
    flint_bitcnt_t bits, min_bits = 2, max_bits = FLINT_BITS;
    nn_ptr pool, res;
    int j;

    if (argc > 1)
        min_bits = FLINT_MAX(2, atoi(argv[1]));
    if (argc > 2)
        max_bits = FLINT_MIN(FLINT_BITS, atoi(argv[2]));

    FLINT_TEST_INIT(state);

    /* warm-up (the first timings of a run can be off) */
    {
        nmod_t mod;
        ulong one = 1, r1;
        slong v = 0;
        nmod_init(&mod, 3);
        time_invert(_nmod_vec_invert_naive, &r1, &one, 1, 1, mod, &v);
    }

    flint_printf("naive / shoup | naive / generic (time ratios)\n");
    flint_printf("bits");
    for (j = 0; j < nlens; j++)
        flint_printf(" %11wu", lens[j]);
    flint_printf("\n");

    pool = _nmod_vec_init(POOL_ENTRIES);
    res = _nmod_vec_init(lens[nlens - 1]);

    for (bits = min_bits; bits <= max_bits; bits++)
    {
        nmod_t mod;
        slong k, v = 0;

        nmod_init(&mod, n_randprime(state, bits, 1));
        for (k = 0; k < POOL_ENTRIES; k++)
            pool[k] = 1 + n_randint(state, mod.n - 1);

        flint_printf("%4wu", bits);
        for (j = 0; j < nlens; j++)
        {
            ulong len = lens[j];
            slong nvec = POOL_ENTRIES / len;
            double tn, ts, tg;

            tn = time_invert(_nmod_vec_invert_naive, res, pool, nvec, len, mod, &v);
            tg = time_invert(_nmod_vec_invert_generic, res, pool, nvec, len, mod, &v);

            if (NMOD_CAN_USE_SHOUP(mod))
            {
                ts = time_invert(_nmod_vec_invert_shoup, res, pool, nvec, len, mod, &v);
                flint_printf(" %5.2f|%5.2f", tn / ts, tn / tg);
            }
            else
                flint_printf("    na|%5.2f", tn / tg);

            v = 0;
        }
        flint_printf("\n");
    }

    _nmod_vec_clear(pool);
    _nmod_vec_clear(res);

    FLINT_TEST_CLEAR(state);
    return 0;
}
