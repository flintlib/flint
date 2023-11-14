/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"
#include "profiler.h"

/* approx nbits(norm(x)) = n */
static void fmpzi_randbits_norm(fmpzi_t x, flint_rand_t state, flint_bitcnt_t n)
{
    fmpz_randbits(fmpzi_realref(x), state, (n+1)/2);
    fmpz_randtest(fmpzi_imagref(x), state, (n+0)/2);
    if (n_randint(state, 2))
        fmpz_swap(fmpzi_realref(x), fmpzi_imagref(x));
}

double profile_it(ulong a_bits, ulong b_bits, ulong g_bits, flint_rand_t state)
{
    fmpzi_t g, a, b;
    slong i, nreps;
    timeit_t timer;

    fmpzi_init(g);
    fmpzi_init(a);
    fmpzi_init(b);

    nreps = 1 + 5000000/(1 + a_bits + b_bits);

    timeit_start(timer);
    for (i = 0; i < nreps; i++)
    {
        fmpzi_randbits_norm(a, state, FLINT_MAX(1, a_bits - g_bits));
        fmpzi_randbits_norm(b, state, FLINT_MAX(1, b_bits - g_bits));
        fmpzi_randbits_norm(g, state, g_bits);
        fmpzi_mul(a, a, g);
        fmpzi_mul(b, b, g);
        fmpzi_gcd_shortest(g, a, b);
        fmpzi_gcd_shortest(g, b, a);
    }
    timeit_stop(timer);

    fmpzi_clear(g);
    fmpzi_clear(a);
    fmpzi_clear(b);

    return timer->wall/nreps;
}


int main(void)
{
    ulong i;
    ulong a_bits, b_bits, g_bits;
    flint_rand_t state;

    flint_randinit(state);

    for (a_bits = 100000; a_bits < 10000000; a_bits += 1 + a_bits/2)
    for (b_bits = a_bits; b_bits < 10000000; b_bits += 1 + b_bits/2)
    {
        flint_printf("a_bits:%8wu, b_bits:%8wu |", a_bits, b_bits);
        for (g_bits = 1; g_bits < a_bits/2; g_bits += 100 + g_bits/2)
        {
            double t = profile_it(a_bits, b_bits, g_bits, state);
            flint_printf(" %5.0f", t);
            fflush(stdout);
        }
        flint_printf("\n");
    }

    flint_randclear(state);

    return 0;
}
