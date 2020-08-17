/*
    Copyright (C) 2012 William Hart
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "mpn_extras.h"

int
fmpz_is_probabprime(const fmpz_t n)
{
    if (!COEFF_IS_MPZ(*n))
    {
        slong v = *n;

        if (v <= 1)
            return 0;
        return n_is_probabprime(v);
    }
    else
    {
        __mpz_struct * z;
        mp_ptr d;
        slong size, bits, trial_primes;

        z = COEFF_TO_PTR(*n);
        size = z->_mp_size;
        d = z->_mp_d;

        if (size < 0)
            return 0;
        if (size == 1)
            return n_is_probabprime(d[0]);

        if (d[0] % 2 == 0)
            return 0;

        bits = size * FLINT_BITS + FLINT_BIT_COUNT(d[size-1]);
        trial_primes = bits;

        if (flint_mpn_factor_trial(d, size, 1, trial_primes))
            return 0;

        if (fmpz_is_square(n))
            return 0;

        return fmpz_is_probabprime_BPSW(n);
    }
}
