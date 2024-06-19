/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

const ulong * n_primes_arr_readonly(ulong num_primes)
{
    slong m;

    if (num_primes < 1)
        return NULL;

    m = FLINT_CLOG2(num_primes);
    if (m >= _flint_primes_used)
        n_compute_primes(num_primes);

    return _flint_primes[m];
}
