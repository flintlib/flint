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
n_primes_clear(n_primes_t iter)
{
    if (iter->small_primes != flint_primes_small)
        flint_free(iter->small_primes);

    if (iter->sieve != NULL)
        flint_free(iter->sieve);
}
