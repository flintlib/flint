/*
    Copyright (C) 2009 Thomas Boothby
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
    
int n_is_oddprime_binary(mp_limb_t n) 
{
    ulong diff, prime_lo, prime_hi;
    const mp_limb_t * primes;

    n_prime_pi_bounds(&prime_lo, &prime_hi, n);
    primes = n_primes_arr_readonly(prime_hi + 1);

    prime_hi--; /* convert to indices of primes in table */
    prime_lo--;

    if (n == primes[prime_hi]) return 1;
    if (n > primes[prime_hi]) return 0;

    diff = (prime_hi - prime_lo + 1) / 2;

    while (1)
    {
        ulong diff2;
        if (primes[prime_lo + diff] <= n) prime_lo += diff;
        if (diff <= UWORD(1)) break;
        diff  = (diff + 1)/2;
        diff2 = (prime_hi - prime_lo + 1)/2;
        if (diff > diff2) diff = diff2;
    }

    return (n == primes[prime_lo]);
}
