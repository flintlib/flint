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
n_primes_init(n_primes_t iter)
{
    iter->small_i = 0;
    iter->small_primes = (unsigned int *) flint_primes_small;
    iter->small_num = FLINT_NUM_PRIMES_SMALL;

    iter->sieve_i = 0;
    iter->sieve_num = 0;
    iter->sieve_a = 0;
    iter->sieve_b = 0;
    iter->sieve = NULL;
}
