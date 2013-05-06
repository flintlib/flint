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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

mp_bitcnt_t fmpz_val2(const fmpz_t x)
{
    fmpz c = *x;
    mp_bitcnt_t t;

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 0)
            t = 0;
        else
            count_trailing_zeros(t, FLINT_ABS(c));
    }
    else
    {
        mp_limb_t *d = (COEFF_TO_PTR(c))->_mp_d;
        mp_bitcnt_t u;

        t = 0;
        while (*d == 0)
        {
            d++;
            t += FLINT_BITS;
        }

        count_trailing_zeros(u, *d);
        t += u;
    }

    return t;
}
