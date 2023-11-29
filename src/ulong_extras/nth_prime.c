/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

mp_limb_t n_nth_prime(ulong n)
{
    if (n == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (n_nth_prime). n_nth_prime(0) is undefined.\n");
    }

    return n_primes_arr_readonly(n)[n-1];
}

void n_nth_prime_bounds(mp_limb_t *lo, mp_limb_t *hi, ulong n)
{
    int bits, ll;
    double llo, lhi;

    /* Lower and upper bounds for ln(n) */
    bits = FLINT_BIT_COUNT(n);
    llo = (bits-1) * 0.6931471;
    lhi = bits * 0.6931472;

    /* Lower bound for ln(ln(n)) */
    if      (n < 16)        ll = 0;
    else if (n < 1619)      ll = 1;
    else if (n < 528491312) ll = 2;
    else                    ll = 3;

    *lo = (mp_limb_t) (n * (llo + ll - 1));
    *hi = (mp_limb_t) (n * (lhi + (ll+1) - (n >= 15985 ? 0.9427 : 0.0)));
}
