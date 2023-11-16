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

void
mag_log1p(mag_t z, const mag_t x)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            mag_zero(z);
        else
            mag_inf(z);
    }
    else
    {
        double t;
        fmpz exp = MAG_EXP(x);

        if (!COEFF_IS_MPZ(exp))
        {
            /* Quick bound by x */
            if (exp < -10)
            {
                mag_set(z, x);
            }
            else if (exp < 1000)
            {
                t = ldexp(MAG_MAN(x), exp - MAG_BITS);
                t = (1.0 + t) * (1 + 1e-14);
                t = mag_d_log_upper_bound(t);
                mag_set_d(z, t);
            }
            else
            {
                /* log(2^(exp-1) * (2*v)) = exp*log(2) + log(2*v), 1 <= 2v < 2 */
                t = (MAG_MAN(x) + 1) * ldexp(1.0, 1 - MAG_BITS);
                t = mag_d_log_upper_bound(t);
                t += (exp - 1) * 0.69314718055994530942;
                t *= (1 + 1e-13);
                mag_set_d(z, t);
            }
        }
        else if (fmpz_sgn(MAG_EXPREF(x)) < 0)
        {
            /* Quick bound by x */
            mag_set(z, x);
        }
        else
        {
            mag_add_ui(z, x, 1);
            mag_log(z, z);
        }
    }
}

