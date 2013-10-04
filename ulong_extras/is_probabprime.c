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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

/*
    This function is used by n_is_prime up to 2^64 and *must* therefore
    act as a primality proof up to that limit. 

    Currently it acts as such all the way up to 2^64.
*/
int n_is_probabprime(mp_limb_t n)
{
    mp_limb_t d;
    unsigned int norm;
	mp_limb_t ninv;

    if (n <= UWORD(1)) return 0;
    if (n == UWORD(2)) return 1;
    if ((n & UWORD(1)) == 0) return 0;

    if (n_is_perfect_power235(n)) return 0;

#if FLINT64
    if (n >= UWORD(10000000000000000)) return n_is_probabprime_BPSW(n);
#endif

    d = n - 1;
    count_trailing_zeros(norm, d);
    d >>= norm;

#if FLINT64
    if (n < UWORD(1122004669633))
#else
    if (n < UWORD(2147483648))
#endif  
    {
        double npre;
        if (n < FLINT_ODDPRIME_SMALL_CUTOFF)
            return n_is_oddprime_small(n);

        if (n < FLINT_PRIMES_TAB_DEFAULT_CUTOFF)
            return n_is_oddprime_binary(n);
      
        npre = n_precompute_inverse(n);

        if (n < UWORD(9080191)) 
        {
            if (n_is_strong_probabprime_precomp(n, npre, UWORD(31), d)
                && n_is_strong_probabprime_precomp(n, npre, UWORD(73), d)) return 1;
            else return 0;
        }

#if FLINT64
        if (n < UWORD(4759123141))
        {
#endif
        if (n_is_strong_probabprime_precomp(n, npre, UWORD(2), d) 
            && n_is_strong_probabprime_precomp(n, npre, UWORD(7), d) 
            && n_is_strong_probabprime_precomp(n, npre, UWORD(61), d)) return 1;
        else return 0;
#if FLINT64
        }

        if (n_is_strong_probabprime_precomp(n, npre, UWORD(2), d) 
            && n_is_strong_probabprime_precomp(n, npre, UWORD(13), d) 
            && n_is_strong_probabprime_precomp(n, npre, UWORD(23), d) 
            && n_is_strong_probabprime_precomp(n, npre, UWORD(1662803), d))
            if (n != UWORD(46856248255981)) return 1;
        return 0;
#endif
    }

	ninv = n_preinvert_limb(n);

    if (n_is_strong_probabprime2_preinv(n, ninv, UWORD(2), d) 
        && n_is_strong_probabprime2_preinv(n, ninv, UWORD(3), d) 
        && n_is_strong_probabprime2_preinv(n, ninv, UWORD(7), d) 
        && n_is_strong_probabprime2_preinv(n, ninv, UWORD(61), d) 
        && n_is_strong_probabprime2_preinv(n, ninv, UWORD(24251), d))
#if FLINT64
        if (n != UWORD(46856248255981)) 
#endif
        return 1;

    return 0;
}
