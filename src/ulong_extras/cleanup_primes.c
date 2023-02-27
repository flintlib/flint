/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

void
n_cleanup_primes()
{
    int i;

    for (i = 0; i < _flint_primes_used; i++)
    {
        if (i < _flint_primes_used - 1 && _flint_primes[i] == _flint_primes[i+1])
            continue;

        flint_free(_flint_primes[i]);
        flint_free(_flint_prime_inverses[i]);
    }

    _flint_primes_used = 0;
}

