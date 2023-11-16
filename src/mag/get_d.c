/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "mag.h"

double
mag_get_d(const mag_t z)
{
    if (mag_is_zero(z))
    {
        return 0.0;
    }
    else if (mag_is_inf(z))
    {
        return D_INF;
    }
    else if (MAG_EXP(z) < -1000 || MAG_EXP(z) > 1000)
    {
        if (fmpz_sgn(MAG_EXPREF(z)) < 0)
            return ldexp(1.0, -1000);
        else
            return D_INF;
    }
    else
    {
        return ldexp(MAG_MAN(z), MAG_EXP(z) - MAG_BITS);
    }
}

