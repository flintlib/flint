/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong n_factor_trial(n_factor_t * factors, ulong n, ulong num_primes)
{
    return n_factor_trial_range(factors, n, UWORD(0), num_primes);
}
