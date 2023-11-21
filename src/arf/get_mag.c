/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_get_mag(mag_t y, const arf_t x)
{
    if (arf_is_zero(x))
    {
        mag_zero(y);
    }
    else if (arf_is_special(x))
    {
        mag_inf(y);
    }
    else
    {
        mp_limb_t t, u;

        ARF_GET_TOP_LIMB(t, x);
        t = (t >> (FLINT_BITS - MAG_BITS)) + LIMB_ONE;

        /* may have rounded up to next power of two */
        u = t >> MAG_BITS;
        /* todo: avoid the addition? check agreement with mag_fast_init_set_arf */
        t = (t >> u) + (u & t);

        _fmpz_add_fast(MAG_EXPREF(y), ARF_EXPREF(x), u);
        MAG_MAN(y) = t;
    }
}

void
arf_get_mag_lower(mag_t y, const arf_t x)
{
    if (arf_is_zero(x))
    {
        mag_zero(y);
    }
    else if (arf_is_special(x))
    {
        if (arf_is_nan(x))
            mag_zero(y);
        else
            mag_inf(y);
    }
    else
    {
        mp_limb_t t;
        ARF_GET_TOP_LIMB(t, x);
        MAG_MAN(y) = t >> (FLINT_BITS - MAG_BITS);
        _fmpz_set_fast(MAG_EXPREF(y), ARF_EXPREF(x));
    }
}

