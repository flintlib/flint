/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"

double
mag_get_d_log2_approx(const mag_t x)
{
    if (mag_is_zero(x))
    {
        return COEFF_MIN;
    }
    else if (mag_is_inf(x))
    {
        return COEFF_MAX;
    }
    else if (COEFF_IS_MPZ(MAG_EXP(x)))
    {
        if (fmpz_sgn(MAG_EXPREF(x)) < 0)
            return COEFF_MIN;
        else
            return COEFF_MAX;
    }
    else
    {
        slong e = MAG_EXP(x);

        if (e < -20 || e > 20)
            return e;
        else
            return e + 1.4426950408889634074 *
                mag_d_log_upper_bound(MAG_MAN(x) * ldexp(1.0, -MAG_BITS));
    }
}

