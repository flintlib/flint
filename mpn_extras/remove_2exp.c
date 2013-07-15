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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"


mp_size_t flint_mpn_remove_2exp(mp_ptr x, mp_size_t xsize, mp_bitcnt_t *bits)
{
    mp_size_t shift_limbs, reduced_size;
    mp_bitcnt_t shift_bits;

    *bits = mpn_scan1(x, 0);

    if (*bits == 0)
        return xsize;

    shift_limbs = *bits / FLINT_BITS;
    shift_bits = *bits % FLINT_BITS;
    reduced_size = xsize - shift_limbs;

    if (shift_bits)
    {
        mpn_rshift(x, x + shift_limbs, reduced_size, shift_bits);
        if (x[reduced_size - 1] == 0)
            reduced_size -= 1;
    }
    else
    {
        flint_mpn_copyi(x, x + shift_limbs, reduced_size);
    }
    return reduced_size;
}
