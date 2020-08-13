/*
    Copyright (C) 2011 Fredrik Johansson

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

void
_fmpz_vec_lcm(fmpz_t res, const fmpz * vec, slong len)
{
    slong i;

    fmpz_one(res);
    for (i = 0; i < len && !fmpz_is_zero(res); i++)
        fmpz_lcm(res, res, vec + i);

    fmpz_abs(res, res);
}
