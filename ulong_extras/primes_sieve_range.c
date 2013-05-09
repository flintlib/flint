/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"

static void
mark(char * sieve, mp_limb_t a, len_t len, mp_limb_t p)
{
    mp_limb_t t;

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

void
n_sieve_odd(char * sieve, ulong n, mp_limb_t a,
    unsigned int * sieve_primes, mp_limb_t bound)
{
    len_t i;
    mp_limb_t p;

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
n_primes_sieve_range(n_primes_t iter, mp_limb_t a, mp_limb_t b)
{
    mp_limb_t bound;
    ulong len, odd_len;

    /* a and b must be odd */
    a += (a % 2 == 0);
    b -= (b % 2 == 0);

    len = b - a + 2;
    odd_len = len / 2;

    if (a < 3 || b < a || len > FLINT_SIEVE_SIZE)
    {
        printf("invalid sieve range %lu,%lu!\n", a, b);
        abort();
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
