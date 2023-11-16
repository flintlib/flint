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
mag_addmul(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_zero(z))
    {
        mag_mul(z, x, y);
    }
    else if (mag_is_inf(z) || mag_is_inf(x) || mag_is_inf(y))
    {
        mag_inf(z);
    }
    else if (mag_is_zero(x) || mag_is_zero(y))
    {
        return;
    }
    else
    {
        slong shift;
        fmpz_t e;

        fmpz_init(e);

        /* x*y < 2^e */
        _fmpz_add2_fast(e, MAG_EXPREF(x), MAG_EXPREF(y), 0);
        shift = _fmpz_sub_small(MAG_EXPREF(z), e);

        if (shift >= 0)
        {
            if (shift >= MAG_BITS)
                MAG_MAN(z)++;
            else
                MAG_MAN(z) = MAG_MAN(z) + (MAG_FIXMUL(MAG_MAN(x),
                    MAG_MAN(y)) >> shift) + LIMB_ONE;
        }
        else
        {
            shift = -shift;
            fmpz_swap(MAG_EXPREF(z), e);

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y))
                    + (2 * LIMB_ONE);
            else
                MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y))
                    + (MAG_MAN(z) >> shift) +  (2 * LIMB_ONE);

            MAG_ADJUST_ONE_TOO_SMALL(z);
        }

        MAG_ADJUST_ONE_TOO_LARGE(z);
        fmpz_clear(e);
    }
}

