/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong
n_primes_next(n_primes_t iter)
{
    if (iter->small_i < iter->small_num)
        return iter->small_primes[(iter->small_i)++];

    for (;;)
    {
        while (iter->sieve_i < iter->sieve_num)
            if (iter->sieve[iter->sieve_i++] != 0)
                return iter->sieve_a + 2 * (iter->sieve_i - 1);

        if (iter->sieve_b == 0)
            n_primes_jump_after(iter, iter->small_primes[iter->small_num-1]);
        else
            n_primes_jump_after(iter, iter->sieve_b);
    }
}
