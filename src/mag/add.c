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
mag_add(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_zero(x))
    {
        mag_set(z, y);
    }
    else if (mag_is_zero(y))
    {
        mag_set(z, x);
    }
    else if (mag_is_inf(x) || mag_is_inf(y))
    {
        mag_inf(z);
    }
    else
    {
        slong shift;
        shift = _fmpz_sub_small(MAG_EXPREF(x), MAG_EXPREF(y));

        if (shift == 0)
        {
            _fmpz_set_fast(MAG_EXPREF(z), MAG_EXPREF(x));
            MAG_MAN(z) = MAG_MAN(x) + MAG_MAN(y);
            MAG_ADJUST_ONE_TOO_LARGE(z); /* may need two adjustments */
        }
        else if (shift > 0)
        {
            _fmpz_set_fast(MAG_EXPREF(z), MAG_EXPREF(x));

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_MAN(x) + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_MAN(x) + (MAG_MAN(y) >> shift) + LIMB_ONE;
        }
        else
        {
            shift = -shift;
            _fmpz_set_fast(MAG_EXPREF(z), MAG_EXPREF(y));

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_MAN(y) + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_MAN(y) + (MAG_MAN(x) >> shift) + LIMB_ONE;
        }

        /* TODO: combine */
        MAG_ADJUST_ONE_TOO_LARGE(z);
    }
}

void
mag_add_lower(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_zero(x))
    {
        mag_set(z, y);
    }
    else if (mag_is_zero(y))
    {
        mag_set(z, x);
    }
    else if (mag_is_inf(x) || mag_is_inf(y))
    {
        mag_inf(z);
    }
    else
    {
        slong shift;
        shift = _fmpz_sub_small(MAG_EXPREF(x), MAG_EXPREF(y));

        if (shift == 0)
        {
            _fmpz_set_fast(MAG_EXPREF(z), MAG_EXPREF(x));
            MAG_MAN(z) = MAG_MAN(x) + MAG_MAN(y);
        }
        else if (shift > 0)
        {
            _fmpz_set_fast(MAG_EXPREF(z), MAG_EXPREF(x));

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_MAN(x);
            else
                MAG_MAN(z) = MAG_MAN(x) + (MAG_MAN(y) >> shift);
        }
        else
        {
            shift = -shift;
            _fmpz_set_fast(MAG_EXPREF(z), MAG_EXPREF(y));

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_MAN(y);
            else
                MAG_MAN(z) = MAG_MAN(y) + (MAG_MAN(x) >> shift);
        }

        if (MAG_MAN(z) >> MAG_BITS)
        {
            MAG_MAN(z) >>= 1;
            _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(z), 1);
        }
    }
}
