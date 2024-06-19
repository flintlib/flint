/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

static void
mark(char * sieve, ulong a, ulong len, ulong p)
{
    ulong t;

    t = p * p;
    if (t >= a)
    {
        t = (t - a) / 2;
    }
    else
    {
        t = p - ((a - p) / 2) % p;
        if (t == p)
            t = 0;
    }

    while (t < len)
    {
        sieve[t] = 0;
        t += p;
    }
}

static void
n_sieve_odd(char * sieve, ulong n, ulong a,
    unsigned int * sieve_primes, ulong bound)
{
    ulong i;
    ulong p;

    for (i = 0; i < n / 2; i++)
        sieve[i] = 1;

    i = 0;
    while (1)
    {
        i++;
        p = sieve_primes[i];
        if (p > bound)
            break;
        mark(sieve, a, n / 2, p);
    }
}

void
n_primes_sieve_range(n_primes_t iter, ulong a, ulong b)
{
    ulong bound;
    ulong len, odd_len;

    /* a and b must be odd */
    a += (a % 2 == 0);
    b -= (b % 2 == 0);

    len = b - a + 2;
    odd_len = len / 2;

    if (a < 3 || b < a || len > FLINT_SIEVE_SIZE)
    {
        flint_throw(FLINT_ERROR, "invalid sieve range %wu,%wu!\n", a, b);
    }

    bound = n_sqrt(b) + 1;

    if (iter->sieve == NULL)
        iter->sieve = flint_malloc(FLINT_SIEVE_SIZE / 2 * sizeof(char));

    n_primes_extend_small(iter, bound);
    n_sieve_odd(iter->sieve, len, a, iter->small_primes, bound);

    iter->sieve_i = 0;
    iter->sieve_num = odd_len;
    iter->sieve_a = a;
    iter->sieve_b = b;
}
