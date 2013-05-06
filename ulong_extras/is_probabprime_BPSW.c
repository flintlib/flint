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

    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_is_probabprime_BPSW(mp_limb_t n)
{
    if (n <= 1UL)
        return 0;

    if ((n & 1UL) == 0UL)
    {
        if (n == 2UL)
            return 1;
        return 0;
    }

    if (((n % 10) == 3) || ((n % 10) == 7))
    {
        if (n_is_probabprime_fermat(n, 2) == 0)
            return 0;

        return n_is_probabprime_fibonacci(n);
    }
    else
    {
        mp_limb_t d;

        d = n - 1UL;
        while ((d & 1UL) == 0UL)
            d >>= 1;

        if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
        {
            double npre = n_precompute_inverse(n);
            if (n_is_strong_probabprime_precomp(n, npre, 2L, d) == 0)
                return 0;
        }
        else
        {
            mp_limb_t ninv = n_preinvert_limb(n);
            if (n_is_strong_probabprime2_preinv(n, ninv, 2L, d) == 0)
                return 0;
        }

        return (n_is_probabprime_lucas(n) == 1);
    }
}
