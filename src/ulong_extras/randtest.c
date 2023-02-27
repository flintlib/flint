/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <limits.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

mp_limb_t n_randtest_bits(flint_rand_t state, int bits)
{
    mp_limb_t m;
    mp_limb_t n;

    m = n_randlimb(state);

    if (m & UWORD(7))
    {
        n = n_randbits(state, bits);
    }
    else
    {
        m >>= 3;

        switch (m & UWORD(7))
        {
            case 0:  n = 0;         break;
            case 1:  n = 1;         break;
            case 2:  n = COEFF_MAX; break;
            case 3:  n = WORD_MAX;  break;
            case 4:  n = UWORD_MAX; break;
            case 5:  n =  (UWORD(1)<<n_randint(state, FLINT_BITS)) 
                        - (UWORD(1)<<n_randint(state, FLINT_BITS));
                                    break;
            case 6:  n =  (UWORD(1)<<n_randint(state, FLINT_BITS));
                                    break;
            case 7:  n = -(UWORD(1)<<n_randint(state, FLINT_BITS));
                                    break;
            default: n = 0;
        }

        if (bits < FLINT_BITS)
           n &= ((UWORD(1)<<bits) - UWORD(1)); /* mask it off */

        if (bits) /* set most significant bit */
           n |= (UWORD(1)<<(bits - 1));
        else
           n = 0;
    }

    return n;
}

mp_limb_t n_randtest(flint_rand_t state)
{
    return n_randtest_bits(state, n_randint(state, FLINT_BITS + 1));
}

mp_limb_t n_randtest_not_zero(flint_rand_t state)
{
    mp_limb_t n;

    while ((n = n_randtest(state)) == 0) ;
    return n;
}
