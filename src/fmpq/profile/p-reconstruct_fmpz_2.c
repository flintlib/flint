/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"
#include "profiler.h"


double profile_it(
    flint_rand_t state,
    slong reps,
    flint_bitcnt_t Nbits,
    flint_bitcnt_t Dbits,
    flint_bitcnt_t mbits,
    int print_it)
{
    slong i, j, k;
    fmpq_t x, y;
    fmpz_t N, D, m, r;
    timeit_t timer;
    double new_time, old_time;
    slong outer_reps = 9;

    fmpq_init(x);
    fmpq_init(y);
    fmpz_init(N);
    fmpz_init(D);
    fmpz_init(m);
    fmpz_init(r);

    new_time = 0;
    old_time = 0;

    k = n_randint(state, 2);

    for (j = 0; j < outer_reps; j++)
    {
        fmpq_zero(y);

        fmpz_randbits(N, state, Nbits);
        fmpz_randbits(D, state, Dbits);
        fmpz_randbits(m, state, mbits);

        fmpz_abs(N, N);
        fmpz_abs(D, D);
        fmpz_abs(m, m);

        fmpz_mul_ui(r, N, 2);

        fmpz_randm(fmpq_numref(x), state, r);
        fmpz_sub(fmpq_numref(x), fmpq_numref(x), N);
        fmpz_randm(fmpq_denref(x), state, D);
        fmpz_add_ui(fmpq_denref(x), fmpq_denref(x), 1);
        fmpq_canonicalise(x);

        fmpz_addmul(m, r, D);
        do fmpz_add_ui(m, m, 1);
        while (!fmpq_mod_fmpz(r, x, m));

        if ((j + k) & 1)
        {
            timeit_start(timer);
            for (i = 0; i < reps; i++)
                _fmpq_reconstruct_fmpz_2(fmpq_numref(y), fmpq_denref(y), r, m, N, D);
            timeit_stop(timer);
            new_time += timer->wall;
        }
        else
        {
            timeit_start(timer);
            for (i = 0; i < reps; i++)
                _fmpq_reconstruct_fmpz_2_naive(fmpq_numref(y), fmpq_denref(y), r, m, N, D);
            timeit_stop(timer);
            old_time += timer->wall;
        }

        if (!fmpq_equal(y, x))
        {
            flint_printf("problem %wd\n", (j + k) & 1);
            printf("x: "); fmpq_print(x); printf("\n");
            printf("y: "); fmpq_print(y); printf("\n");
            flint_abort();
        }
    }

    new_time = new_time / reps / outer_reps;
    old_time = old_time / reps / outer_reps;

    if (print_it)
    flint_printf("bits (%wu, %wu, %wu): new %f ms , old %f ms | old/new: %f\n",
                   Nbits, Dbits, mbits, new_time, old_time, old_time/new_time);

    fmpq_clear(x);
    fmpq_clear(y);
    fmpz_clear(N);
    fmpz_clear(D);
    fmpz_clear(m);
    fmpz_clear(r);

    return old_time/new_time;
}

int
main(void)
{
    double total;
    slong i, j, count;
    flint_rand_t state;

    flint_printf("\n");
    fflush(stdout);

    flint_randinit(state);

printf("---- balanced ----\n");

    for (j = 1; j <= 10; j++)
    {
        total = 0;
        count = 0;
        i = (j < 2) ? 10 : FLINT_BITS*(j - 1) + 1;
        for (; i <= FLINT_BITS*j; i += 1)
        {
            total += profile_it(state, 200000/i, i/2, i/2, i, 0);
            count++;
        }
        flint_printf("size %wd average speedup ratio: %f\n", j, total/count);
    }

    profile_it(state,  100,    500,    500,   1000,1);
    profile_it(state,  100,    500,    500,   1000,1);
    profile_it(state,   80,   1000,   1000,   2000,1);
    profile_it(state,   80,   1000,   1000,   2000,1);
    profile_it(state,   50,   2000,   2000,   4000,1);
    profile_it(state,   50,   2000,   2000,   4000,1);
    profile_it(state,   10,  10000,  10000,  20000,1);
    profile_it(state,   10,  10000,  10000,  20000,1);

    profile_it(state,    8,  20000,  20000,  40000,1);
    profile_it(state,    8,  20000,  20000,  40000,1);
    profile_it(state,    6,  40000,  40000,  80000,1);
    profile_it(state,    6,  40000,  40000,  80000,1);
    profile_it(state,    4,  80000,  80000, 160000,1);
    profile_it(state,    1, 160000, 160000, 320000,1);

printf("---- imbalanced ----\n");

    profile_it(state, 500,  50, 30, 80,1);
    profile_it(state, 500,  30, 50, 80,1);

    profile_it(state, 200,  50,100,150,1);
    profile_it(state, 200, 100, 50,150,1);

    profile_it(state, 100, 100,200,300,1);
    profile_it(state, 100, 200,100,300,1);

    profile_it(state, 100, 200,400,600,1);
    profile_it(state, 100, 400,200,600,1);

    profile_it(state, 100, 400,800,1200,1);
    profile_it(state, 100, 800,400,1200,1);

    profile_it(state, 100, 1000,2000,3000,1);
    profile_it(state, 100, 2000,1000,3000,1);

    profile_it(state,  80, 2000,4000,6000,1);
    profile_it(state,  80, 4000,2000,6000,1);

    profile_it(state,  30, 4000,8000,12000,1);
    profile_it(state,  30, 8000,4000,12000,1);

    profile_it(state,  10, 10000,20000,30000,1);
    profile_it(state,  10, 20000,10000,30000,1);

    profile_it(state,   4, 20000,40000,60000,1);
    profile_it(state,   4, 40000,20000,60000,1);


    flint_randclear(state);

    flint_cleanup_master();
    return 0;
}

