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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <limits.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

mp_limb_signed_t z_randtest(flint_rand_t state)
{
    mp_limb_t m;
    mp_limb_signed_t z;

    m = n_randlimb(state);

    if (m & 7UL)
    {
        z = n_randbits(state, n_randint(state, FLINT_BITS));
    }
    else
    {
        m >>= 3;

        switch (m % 7UL)
        {
            case 0:  z = 0;         break;
            case 1:  z = 1;         break;
            case 2:  z = -1;        break;
            case 3:  z = COEFF_MAX; break;
            case 4:  z = COEFF_MIN; break;
            case 5:  z = LONG_MAX;  break;
            case 6:  z = LONG_MIN;  break;
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

