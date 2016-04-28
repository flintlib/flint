/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void
fmpz_mods(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c2 = *h;

    if (!COEFF_IS_MPZ(c2))      /* h is small */
    {
        ulong tmp;

        tmp = FLINT_ABS(c2);

        fmpz_mod(f, g, h);
        if (fmpz_cmp_ui(f, tmp / 2) > 0)
        {
            fmpz_sub_ui(f, f, tmp);
        }
    }
    else                        /* h is large */
    {
        fmpz_t tmp, rtmp;

        fmpz_init(tmp);
        fmpz_init(rtmp);
        fmpz_abs(tmp, h);
        fmpz_fdiv_q_2exp(rtmp, tmp, 1);

        fmpz_mod(f, g, h);
        if (fmpz_cmp(f, rtmp) > 0)
        {
            fmpz_sub(f, f, tmp);
        }
        fmpz_clear(tmp);
        fmpz_clear(rtmp);
    }
}
