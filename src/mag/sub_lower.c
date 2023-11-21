/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"
#include "arf.h"

void
mag_sub_lower(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_special(x) || mag_is_special(y))
    {
        if (mag_is_zero(y))
            mag_set(z, x);
        else if (mag_is_zero(x) || mag_is_inf(y))
            mag_zero(z);
        else
            mag_inf(z);
    }
    else
    {
        slong shift, fix;

        shift = _fmpz_sub_small(MAG_EXPREF(x), MAG_EXPREF(y));

        /* y > x */
        if (shift < 0)
        {
            mag_zero(z);
        }
        else if (shift == 0)
        {
            if (MAG_MAN(y) >= MAG_MAN(x))
            {
                mag_zero(z);
            }
            else
            {
                MAG_MAN(z) = MAG_MAN(x) - MAG_MAN(y);
                fix = MAG_BITS - FLINT_BIT_COUNT(MAG_MAN(z));
                MAG_MAN(z) <<= fix;
                _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(x), -fix);
            }
        }
        else
        {
            if (shift <= MAG_BITS)
            {
                mp_limb_t c = MAG_MAN(x) - (MAG_MAN(y) >> shift) - 1;

                /* too much cancellation -- compute precisely */
                if (c < (UWORD(1) << (MAG_BITS - 4)))
                {
                    arf_t t, u;
                    arf_init(t);
                    arf_init(u);
                    arf_set_mag(t, x);
                    arf_set_mag(u, y);
                    arf_sub(t, t, u, MAG_BITS, ARF_RND_DOWN);
                    arf_get_mag_lower(z, t);
                    arf_clear(t);
                    arf_clear(u);
                    return;
                }

                MAG_MAN(z) = c;
            }
            else
            {
                MAG_MAN(z) = MAG_MAN(x) - 1;
            }

            fix = MAG_BITS - FLINT_BIT_COUNT(MAG_MAN(z));
            MAG_MAN(z) <<= fix;
            _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(x), -fix);
        }
    }
}

