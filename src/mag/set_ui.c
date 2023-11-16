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
mag_set_ui(mag_t z, ulong x)
{
    _fmpz_demote(MAG_EXPREF(z));

    if (x == 0)
    {
        MAG_EXP(z) = 0;
        MAG_MAN(z) = 0;
    }
    else
    {
        slong bits;
        mp_limb_t overflow;

        bits = flint_clz(x);
        bits = FLINT_BITS - bits;

        if (bits <= MAG_BITS)
        {
            x = x << (MAG_BITS - bits);
        }
        else
        {
            x = (x >> (bits - MAG_BITS)) + 1;
            overflow = x >> MAG_BITS;
            bits += overflow;
            x >>= overflow;
        }

        MAG_EXP(z) = bits;
        MAG_MAN(z) = x;
    }
}

void
mag_set_ui_lower(mag_t z, ulong x)
{
    _fmpz_demote(MAG_EXPREF(z));

    if (x == 0)
    {
        MAG_EXP(z) = 0;
        MAG_MAN(z) = 0;
    }
    else
    {
        unsigned int bits;

        bits = flint_clz(x);
        bits = FLINT_BITS - bits;

        if (bits <= MAG_BITS)
            x <<= (MAG_BITS - bits);
        else
            x >>= (bits - MAG_BITS);

        MAG_EXP(z) = bits;
        MAG_MAN(z) = x;
    }
}

