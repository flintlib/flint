/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

int
mag_cmp_2exp_si(const mag_t x, slong e)
{
    int ispow2;

    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            return -1;
        return 1;
    }

    ispow2 = (MAG_MAN(x) == MAG_ONE_HALF);

    /* Fast path. */
    if (!COEFF_IS_MPZ(MAG_EXP(x)))
    {
        if (ispow2 && (MAG_EXP(x) - 1 == e))
            return 0;
        else
            return (MAG_EXP(x) <= e) ? -1 : 1;
    }


    if (ispow2)
    {
        fmpz_t t;
        fmpz_init(t);

        fmpz_one(t);
        fmpz_add_si(t, t, e);

        if (fmpz_equal(MAG_EXPREF(x), t))
        {
            fmpz_clear(t);
            return 0;
        }

        fmpz_clear(t);
    }

    return (fmpz_cmp_si(MAG_EXPREF(x), e) <= 0) ? -1 : 1;
}

