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
mag_expm1(mag_t y, const mag_t x)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            mag_zero(y);
        else
            mag_inf(y);
    }
    else if (COEFF_IS_MPZ(MAG_EXP(x)))
    {
        if (fmpz_sgn(MAG_EXPREF(x)) > 0)
        {
            mag_inf(y);
        }
        else
        {
            fmpz_set(MAG_EXPREF(y), MAG_EXPREF(x));
            MAG_MAN(y) = MAG_MAN(x) + 1;
            MAG_ADJUST_ONE_TOO_LARGE(y);
        }
    }
    else
    {
        slong e = MAG_EXP(x);

        if (e <= -16)
        {
            fmpz_set(MAG_EXPREF(y), MAG_EXPREF(x));

            if (e < -MAG_BITS)
                MAG_MAN(y) = MAG_MAN(x) + 1;
            else
                MAG_MAN(y) = MAG_MAN(x) + (MAG_ONE_HALF >> 15);

            MAG_ADJUST_ONE_TOO_LARGE(y);
        }
        else
        {
            mag_exp(y, x);

            /* subtract 1 */
            if (e < 6)
            {
                slong e2 = MAG_EXP(y);
                unsigned int c;

                if (e2 < MAG_BITS)
                {
                    /* this is correct since exp(x) >= 1 */
                    MAG_MAN(y) -= (MAG_ONE_HALF >> (e2 - 1));
                    c = MAG_BITS - FLINT_BIT_COUNT(MAG_MAN(y));
                    MAG_MAN(y) <<= c;
                    MAG_EXP(y) -= c;
                }
            }
        }
    }
}

