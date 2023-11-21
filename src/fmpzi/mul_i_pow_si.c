/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

void
fmpzi_mul_i_pow_si(fmpzi_t res, const fmpzi_t z, slong k)
{
    k &= 3;

    if (k == 0)
        fmpzi_set(res, z);
    else if (k == 1)
        fmpzi_mul_i(res, z);
    else if (k == 2)
        fmpzi_neg(res, z);
    else
        fmpzi_div_i(res, z);
}
