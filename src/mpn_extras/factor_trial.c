/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"
#include "ulong_extras.h"


int flint_mpn_factor_trial(mp_srcptr x, mp_size_t xsize, slong start, slong stop)
{
    slong i;
    const mp_limb_t * primes;

    primes = n_primes_arr_readonly(stop);

    for (i = start; i < stop; i++)
    {
        if (flint_mpn_divisible_1_p(x, xsize, primes[i]))
            return i;
    }
    return 0;
}
