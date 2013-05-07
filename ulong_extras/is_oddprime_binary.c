/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009 Thomas Boothby
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
    
int n_is_oddprime_binary(mp_limb_t n) 
{
    ulong diff, prime_lo, prime_hi;

    n_prime_pi_bounds(&prime_lo, &prime_hi, n);
    if (prime_hi >= flint_num_primes) prime_hi = flint_num_primes - 1;

    if (n == flint_primes[prime_hi]) return 1;
    if (n > flint_primes[prime_hi]) return 0;
    
    diff = (prime_hi - prime_lo + 1) / 2;

    while (1)
    {
        ulong diff2;
        if (flint_primes[prime_lo + diff] <= n) prime_lo += diff;
        if (diff <= 1UL) break;
        diff  = (diff + 1)/2;
        diff2 = (prime_hi - prime_lo + 1)/2;
        if (diff > diff2) diff = diff2;
    }
       
    return (n == flint_primes[prime_lo]);
}
