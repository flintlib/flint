/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

/*
    Perfect square test. Small values use n_is_square. Multi-limb values
    go through flint_mpn_is_square: quadratic residue screening (modulo 256
    and the factors of 2^48 - 1) rejects nearly all random nonsquares in
    O(n) limb operations, and the remaining candidates are certified by
    computing the square root (Newton iteration with fft_small
    multiplication for large inputs).
*/
int
fmpz_is_square(const fmpz_t x)
{
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        if (c <= 1)
            return (c >= 0);
        return n_is_square(c);
    }
    else
    {
        mpz_srcptr mx = COEFF_TO_PTR(c);
        mp_size_t an = mx->_mp_size;

        if (an < 0)
            return 0;

        return flint_mpn_is_square(mx->_mp_d, an);
    }
}
