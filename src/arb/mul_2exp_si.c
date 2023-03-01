/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* mag_mul_2exp_si is non-inline, but avoid overhead here */
static __inline__ void
_mag_mul_2exp_si(mag_t z, const mag_t x, slong y)
{
    if (mag_is_special(x))
    {
        mag_set(z, x);
    }
    else
    {
        if (y >= ADD2_FAST_MIN && y <= ADD2_FAST_MAX)
            _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(x), y);
        else
            fmpz_add_si(MAG_EXPREF(z), MAG_EXPREF(x), y);
        MAG_MAN(z) = MAG_MAN(x);
    }
}

void
arb_mul_2exp_si(arb_t y, const arb_t x, slong e)
{
    arf_mul_2exp_si(arb_midref(y), arb_midref(x), e);
    _mag_mul_2exp_si(arb_radref(y), arb_radref(x), e);
}

