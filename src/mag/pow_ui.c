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
mag_pow_ui(mag_t z, const mag_t x, ulong e)
{
    if (mag_is_inf(x))
    {
        mag_inf(z);
    }
    else if (e <= 2)
    {
        if (e == 0)
            mag_one(z);
        else if (e == 1)
            mag_set(z, x);
        else
            mag_mul(z, x, x);
    }
    else
    {
        mag_t y;
        int i, bits;

        mag_init_set(y, x);

        bits = FLINT_BIT_COUNT(e);

        for (i = bits - 2; i >= 0; i--)
        {
            mag_mul(y, y, y);
            if (e & (UWORD(1) << i))
                mag_mul(y, y, x);
        }

        mag_swap(z, y);
        mag_clear(y);
    }
}

void
mag_pow_ui_lower(mag_t z, const mag_t x, ulong e)
{
    if (e <= 2)
    {
        if (e == 0)
            mag_one(z);
        else if (e == 1)
            mag_set(z, x);
        else
            mag_mul_lower(z, x, x);
    }
    else if (mag_is_inf(x))
    {
        mag_inf(z);
    }
    else
    {
        mag_t y;
        int i, bits;

        mag_init_set(y, x);

        bits = FLINT_BIT_COUNT(e);

        for (i = bits - 2; i >= 0; i--)
        {
            mag_mul_lower(y, y, y);
            if (e & (UWORD(1) << i))
                mag_mul_lower(y, y, x);
        }

        mag_swap(z, y);
        mag_clear(y);
    }
}

