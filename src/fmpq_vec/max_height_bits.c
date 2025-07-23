/*
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

flint_bitcnt_t
_fmpq_vec_max_height_bits(const fmpq * x, slong len)
{
    flint_bitcnt_t bits, max_bits = 0;

    for (slong i = 0; i < len; i++)
    {
        bits = fmpq_height_bits(x + i);
        if (bits > max_bits)
            max_bits = bits;
    }

    return max_bits;
}
