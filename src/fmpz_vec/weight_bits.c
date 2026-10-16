/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

slong
_fmpz_vec_weight_bits(const fmpz * vec, slong len)
{
    slong i;
    ulong w = 0;

    for (i = 0; i < len; i++)
        w += fmpz_bits(vec + i);

    return w;
}
