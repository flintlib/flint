/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_div_lower(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_special(x) || mag_is_special(y))
    {
        if (mag_is_inf(x) && !mag_is_inf(y))
            mag_inf(z);
        else
            mag_zero(z);
    }
    else
    {
        mp_limb_t q;
        slong fix;

#if FLINT_BITS == 64
        q = (MAG_MAN(x) << MAG_BITS) / MAG_MAN(y);
#else
        mp_limb_t hi, lo, r;
        lo = MAG_MAN(x) << MAG_BITS;
        hi = MAG_MAN(x) >> (FLINT_BITS - MAG_BITS);
        udiv_qrnnd(q, r, hi, lo, MAG_MAN(y));
#endif
        fix = q >> MAG_BITS;
        q = q >> fix;

        /* could overflow to next power of two */
        if (q >> MAG_BITS)
        {
            q >>= 1;
            fix += 1;
        }

        MAG_MAN(z) = q;
        _fmpz_sub2_fast(MAG_EXPREF(z), MAG_EXPREF(x), MAG_EXPREF(y), fix);
    }
}

