/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Dana Jacobsen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
    int isprime;
#if FLINT64
    double npre;
#else
    mp_limb_t ninv;
#endif

    if (n <= UWORD(1)) return 0;
    if (n == UWORD(2)) return 1;
    if ((n & UWORD(1)) == 0) return 0;

    if (n < FLINT_ODDPRIME_SMALL_CUTOFF)
        return n_is_oddprime_small(n);
    if (n < FLINT_PRIMES_TAB_DEFAULT_CUTOFF)
        return n_is_oddprime_binary(n);

#if FLINT64
    /* Avoid the unnecessary inverse */
    if (n >= UWORD(1050535501))
        return n_is_probabprime_BPSW(n);
#endif

    isprime = 0;
    d = n - 1;
    count_trailing_zeros(norm, d);
    d >>= norm;

#if !FLINT64

    /* For 32-bit, just the 2-base or 3-base Miller-Rabin is enough */
    /* The preinv functions are faster on 32-bit, and work up to
       2^32 (precomp only works up to 2^31) */
    ninv = n_preinvert_limb(n);

    if (n < UWORD(9080191))
    {
        isprime = n_is_strong_probabprime2_preinv(n, ninv, UWORD(31), d)
               && n_is_strong_probabprime2_preinv(n, ninv, UWORD(73), d);
    }
    else
    {
        isprime = n_is_strong_probabprime2_preinv(n, ninv, UWORD(2), d)
               && n_is_strong_probabprime2_preinv(n, ninv, UWORD(7), d)
               && n_is_strong_probabprime2_preinv(n, ninv, UWORD(61), d);
    }
#else
    npre = n_precompute_inverse(n);

    /* For 64-bit, BPSW seems to be a little bit faster than 3 bases. */
    if (n < UWORD(341531))
    {
        isprime = n_is_strong_probabprime_precomp(n, npre, UWORD(9345883071009581737), d);
    }
    else if (n < UWORD(1050535501))
    {
        isprime = n_is_strong_probabprime_precomp(n, npre, UWORD(336781006125), d)
               && n_is_strong_probabprime_precomp(n, npre, UWORD(9639812373923155), d);
    }
#if 0
    else if (n < UWORD(350269456337))
    {
        isprime = n_is_strong_probabprime_precomp(n, npre, UWORD(4230279247111683200), d)
               && n_is_strong_probabprime_precomp(n, npre, UWORD(14694767155120705706), d)
               && n_is_strong_probabprime_precomp(n, npre, UWORD(16641139526367750375), d);
    }
#endif
    else
    {
        isprime = n_is_probabprime_BPSW(n);
    }

#endif

    return isprime;
}
