/*
    Copyright 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_vec.h"
#include "flint/profiler.h"

void sample_ndiv_qr(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t q, r, a, b;
    fmpz_t nmax;

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, FLINT_BITS);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        if (n_randint(state, 2))
            fmpz_neg(a, a);
        if (n_randint(state, 2))
            fmpz_neg(b, b);
        if (fmpz_is_zero(b))
            fmpz_one(b);
        fmpz_ndiv_qr(q, r, a, b);
    }
    prof_stop();

    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

void sample_fdiv_qr(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t q, r, a, b;
    fmpz_t nmax;

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, FLINT_BITS);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        if (n_randint(state, 2))
            fmpz_neg(a, a);
        if (n_randint(state, 2))
            fmpz_neg(b, b);
        if (fmpz_is_zero(b))
            fmpz_one(b);
        fmpz_fdiv_qr(q, r, a, b);
    }
    prof_stop();

    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

void sample_cdiv_qr(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t q, r, a, b;
    fmpz_t nmax;

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, FLINT_BITS);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        if (n_randint(state, 2))
            fmpz_neg(a, a);
        if (n_randint(state, 2))
            fmpz_neg(b, b);
        if (fmpz_is_zero(b))
            fmpz_one(b);
        fmpz_cdiv_qr(q, r, a, b);
    }
    prof_stop();

    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

void sample_tdiv_qr(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t q, r, a, b;
    fmpz_t nmax;

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, FLINT_BITS);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        if (n_randint(state, 2))
            fmpz_neg(a, a);
        if (n_randint(state, 2))
            fmpz_neg(b, b);
        if (fmpz_is_zero(b))
            fmpz_one(b);
        fmpz_tdiv_qr(q, r, a, b);
    }
    prof_stop();

    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

int main(void)
{
    double min, max;

    prof_repeat(&min, &max, sample_ndiv_qr, NULL);
    flint_printf("fmpz_ndiv_qr:\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    prof_repeat(&min, &max, sample_fdiv_qr, NULL);
    flint_printf("fmpz_fdiv_qr:\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    prof_repeat(&min, &max, sample_cdiv_qr, NULL);
    flint_printf("fmpz_cdiv_qr:\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    prof_repeat(&min, &max, sample_tdiv_qr, NULL);
    flint_printf("fmpz_tdiv_qr:\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    {
        flint_bitcnt_t abits, bbits;
        fmpz * as, * bs;
        fmpz_t q, r;
        slong i, len = 2000000;
        flint_rand_t state;
        timeit_t timer;

        flint_randinit(state);

        as = _fmpz_vec_init(len);
        bs = _fmpz_vec_init(len);

        fmpz_init(q);
        fmpz_init(r);

        for (bbits = 10; bbits <= 300; bbits += bbits/4)
        {
            flint_printf("%4wu: ", bbits);
            for (abits = 10; abits <= 300; abits += abits/4)
            {
                for (i = 0; i < len; i++)
                {
                    fmpz_randbits(as + i, state, abits);
                    fmpz_randbits(bs + i, state, bbits);
                }

                timeit_start(timer);
                for (i = 0; i < len; i++)
                    fmpz_ndiv_qr(q, r, as + i, bs + i);
                timeit_stop(timer);

                flint_printf(" %4wd", timer->wall);
                fflush(stdout);
            }
            flint_printf("\n");
        }

        fmpz_clear(q);
        fmpz_clear(r);

        _fmpz_vec_clear(as, len);
        _fmpz_vec_clear(bs, len);

        flint_randclear(state);
    }

    return 0;
}
