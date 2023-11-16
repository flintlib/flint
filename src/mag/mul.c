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
mag_mul(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_special(x) || mag_is_special(y))
    {
        if (mag_is_inf(x) || mag_is_inf(y))
            mag_inf(z);
        else
            mag_zero(z);
    }
    else
    {
        slong fix;

        MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + LIMB_ONE;
        fix = !(MAG_MAN(z) >> (MAG_BITS - 1));
        MAG_MAN(z) <<= fix;
        _fmpz_add2_fast(MAG_EXPREF(z), MAG_EXPREF(x), MAG_EXPREF(y), -fix);
    }
}

void
mag_mul_lower(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_special(x) || mag_is_special(y))
    {
        if (mag_is_zero(x) || mag_is_zero(y))
            mag_zero(z);
        else
            mag_inf(z);
    }
    else
    {
        slong fix;

        MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y));
        fix = !(MAG_MAN(z) >> (MAG_BITS - 1));
        MAG_MAN(z) <<= fix;
        _fmpz_add2_fast(MAG_EXPREF(z), MAG_EXPREF(x), MAG_EXPREF(y), -fix);
    }
}

