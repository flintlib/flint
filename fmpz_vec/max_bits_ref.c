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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

long
_fmpz_vec_max_bits_ref(const fmpz * vec, long len)
{
    long i, bits, max_bits = 0, sign = 1;

    for (i = 0; i < len; i++)
    {
        bits = fmpz_bits(vec + i);
        if (bits > max_bits)
            max_bits = bits;
        if (fmpz_sgn(vec + i) < 0)
            sign = -1L;
    }

    return max_bits * sign;
}
