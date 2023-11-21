/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_mul_2exp_fmpz(mag_t z, const mag_t x, const fmpz_t y)
{
    if (mag_is_special(x))
    {
        mag_set(z, x);
    }
    else
    {
        _fmpz_add2_fast(MAG_EXPREF(z), MAG_EXPREF(x), y, 0);
        MAG_MAN(z) = MAG_MAN(x);
    }
}

