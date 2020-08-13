/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"

void
n_primes_extend_small(n_primes_t iter, mp_limb_t bound)
{
    while (iter->small_primes[iter->small_num - 2] < bound)
    {
        n_primes_t iter2;
        slong i, num;

        num = iter->small_num * 2;

        if (iter->small_primes == flint_primes_small)
            iter->small_primes = flint_malloc(num * sizeof(unsigned int));
        else
            iter->small_primes = flint_realloc(iter->small_primes,
                num * sizeof(unsigned int));

        n_primes_init(iter2);
        for (i = 0; i < num; i++)
            iter->small_primes[i] = n_primes_next(iter2);
        n_primes_clear(iter2);

        iter->small_num = num;
        iter->small_i = num;
    }
}
