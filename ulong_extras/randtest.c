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
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <limits.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

mp_limb_t n_randtest(flint_rand_t state)
{
    mp_limb_t m;
    mp_limb_t n;

    m = n_randlimb(state);

    if (m & 7UL)
    {
        n = n_randbits(state, n_randint(state, FLINT_BITS + 1));
    }
    else
    {
        m >>= 3;

        switch (m % 5UL)
        {
            case 0:  n = 0;         break;
            case 1:  n = 1;         break;
            case 2:  n = COEFF_MAX; break;
            case 3:  n = LONG_MAX;  break;
            case 4:  n = ULONG_MAX; break;
            default: n = 0;
        }
    }

    return n;
}

mp_limb_t n_randtest_not_zero(flint_rand_t state)
{
    mp_limb_t n;

    while ((n = n_randtest(state)) == 0) ;
    return n;
}
