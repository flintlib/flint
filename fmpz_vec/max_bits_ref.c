/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

slong
_fmpz_vec_max_bits_ref(const fmpz * vec, slong len)
{
    slong i, bits, max_bits = 0, sign = 1;

    for (i = 0; i < len; i++)
    {
        bits = fmpz_bits(vec + i);
        if (bits > max_bits)
            max_bits = bits;
        if (fmpz_sgn(vec + i) < 0)
            sign = WORD(-1);
    }

    return max_bits * sign;
}
