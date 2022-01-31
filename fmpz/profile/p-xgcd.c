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
#include "flint/profiler.h"

void sample_xgcd_small(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t d, x, y, a, b;
    fmpz_t nmax;

    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, SMALL_FMPZ_BITCOUNT_MAX);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        fmpz_xgcd(d, x, y, a, b);
    }
    prof_stop();

    fmpz_clear(d);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

void sample_xgcd_mixed(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t d, x, y, a, b;
    fmpz_t nmax;

    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, FLINT_BITS);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        fmpz_xgcd(d, x, y, a, b);
    }
    prof_stop();

    fmpz_clear(d);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

void sample_xgcd_big(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t d, x, y, a, b;
    fmpz_t nmax;

    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, 512);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        fmpz_xgcd(d, x, y, a, b);
    }
    prof_stop();

    fmpz_clear(d);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

void sample_xgcd_canonical_bezout_small(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t d, x, y, a, b;
    fmpz_t nmax;

    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, SMALL_FMPZ_BITCOUNT_MAX);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        fmpz_xgcd_canonical_bezout(d, x, y, a, b);
    }
    prof_stop();

    fmpz_clear(d);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

void sample_xgcd_canonical_bezout_mixed(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t d, x, y, a, b;
    fmpz_t nmax;

    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, FLINT_BITS);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        fmpz_xgcd_canonical_bezout(d, x, y, a, b);
    }
    prof_stop();

    fmpz_clear(d);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

void sample_xgcd_canonical_bezout_big(void * arg, ulong count)
{
    FLINT_TEST_INIT(state);
    fmpz_t d, x, y, a, b;
    fmpz_t nmax;

    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(a);
    fmpz_init(b);
    
    fmpz_init(nmax);
    fmpz_set_d_2exp(nmax, 1.0, 512);

    prof_start();
    for (int ix = 0; ix < count; ix++)
    {
        fmpz_randm(a, state, nmax);
        fmpz_randm(b, state, nmax);
        fmpz_xgcd_canonical_bezout(d, x, y, a, b);
    }
    prof_stop();

    fmpz_clear(d);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(nmax);

    flint_randclear(state);
}

int main(void)
{
    double min, max;

    prof_repeat(&min, &max, sample_xgcd_small, NULL);
    flint_printf("fmpz_xgcd (small size):\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    prof_repeat(&min, &max, sample_xgcd_mixed, NULL);
    flint_printf("fmpz_xgcd (mixed size):\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    prof_repeat(&min, &max, sample_xgcd_big, NULL);
    flint_printf("fmpz_xgcd (big size):\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    prof_repeat(&min, &max, sample_xgcd_canonical_bezout_small, NULL);
    flint_printf("fmpz_xgcd_canonical_bezout (small size):\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    prof_repeat(&min, &max, sample_xgcd_canonical_bezout_mixed, NULL);
    flint_printf("fmpz_xgcd_canonical_bezout (mixed size):\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

    prof_repeat(&min, &max, sample_xgcd_canonical_bezout_big, NULL);
    flint_printf("fmpz_xgcd_canonical_bezout (big size):\n"
                 "  min time is %.3f cycles\n"
                 "  max time is %.3f cycles\n\n",
                 (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100,
                 (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);
    return 0;
}
