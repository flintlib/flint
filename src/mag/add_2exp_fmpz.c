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
mag_add_2exp_fmpz(mag_t z, const mag_t x, const fmpz_t e)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
        {
            MAG_MAN(z) = MAG_ONE_HALF;
            _fmpz_add_fast(MAG_EXPREF(z), e, 1);
        }
        else
        {
            mag_inf(z);
        }
    }
    else
    {
        slong shift;

        shift = _fmpz_sub_small(MAG_EXPREF(x), e);

        if (shift > 0)
        {
            _fmpz_set_fast(MAG_EXPREF(z), MAG_EXPREF(x));

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_MAN(x) + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_MAN(x) + (LIMB_ONE << (MAG_BITS - shift));
        }
        else
        {
            shift = -shift;

            _fmpz_add_fast(MAG_EXPREF(z), e, 1);

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_ONE_HALF + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_ONE_HALF + (MAG_MAN(x) >> (shift + 1)) + LIMB_ONE;
        }

        MAG_ADJUST_ONE_TOO_LARGE(z);
    }
}

