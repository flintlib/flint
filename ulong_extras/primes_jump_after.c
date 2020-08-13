/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "ulong_extras.h"

void
n_primes_jump_after(n_primes_t iter, mp_limb_t n)
{
    if (n < iter->small_primes[iter->small_num - 1])
    {
        iter->small_i = n_prime_pi(n);
        iter->sieve_a = iter->sieve_b = iter->sieve_num = 0;
    }
    else
    {
        iter->small_i = iter->small_num;
        n_primes_sieve_range(iter, n + 1,
            n + 1 + FLINT_MIN(n, FLINT_SIEVE_SIZE) - 2);
    }
}
