/*
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
#include "long_extras.h"

mp_limb_signed_t z_randtest(flint_rand_t state)
{
    mp_limb_t m;
    mp_limb_signed_t z;

    m = n_randlimb(state);

    if (m & UWORD(7))
    {
        z = n_randbits(state, n_randint(state, FLINT_BITS));
    }
    else
    {
        m >>= 3;

        switch (m % UWORD(7))
        {
            case 0:  z = 0;         break;
            case 1:  z = 1;         break;
            case 2:  z = -1;        break;
            case 3:  z = COEFF_MAX; break;
            case 4:  z = COEFF_MIN; break;
            case 5:  z = WORD_MAX;  break;
            case 6:  z = WORD_MIN;  break;
            default: z = 0;
        }
    }

    return z;
}

mp_limb_signed_t z_randtest_not_zero(flint_rand_t state)
{
    mp_limb_signed_t z;

    while ((z = z_randtest(state)) == 0) ;
    return z;
}

