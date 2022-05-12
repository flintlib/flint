/*
    Copyright 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"

/* negation mod 2^(FLINT_BITS*limbs) + 1 assuming normalized input */
void mpn_negmod_2expp1(mp_limb_t* z, const mp_limb_t* a, mp_size_t limbs)
{
    if (a[limbs] != 0)
    {
        FLINT_ASSERT(a[limbs] == 1);
        z[0] = 1;
        flint_mpn_zero(z + 1, limbs);
    }
    else
    {
        mpn_com(z, a, limbs);
        z[limbs] = mpn_add_1(z, z, limbs, 2);
        if (z[limbs] != 0)
        {
            if (z[0] != 0)
            {
                FLINT_ASSERT(z[0] == 1);
                z[0] = 0;
                z[limbs] = 0;
            }
        }
    }
}

