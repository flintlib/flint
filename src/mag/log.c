/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
#else
# include <math.h>
#endif

void
mag_neg_log(mag_t z, const mag_t x)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            mag_inf(z);
        else
            mag_zero(z);
    }
    else
    {
        double t;
        fmpz exp = MAG_EXP(x);

        if (!COEFF_IS_MPZ(exp))
        {
            if (exp >= 1)
            {
                mag_zero(z);
            }
            else if (exp > -(1000 - MAG_BITS))
            {
                t = ldexp(MAG_MAN(x), exp - MAG_BITS);
                t = -mag_d_log_lower_bound(t);
                mag_set_d(z, t);
            }
            else
            {
                /* -log(2^(exp-1) * (2*v)) = -exp*log(2) - log(2*v), 1 <= 2v < 2 */
                t = MAG_MAN(x) * ldexp(1.0, 1 - MAG_BITS);
                t = -mag_d_log_lower_bound(t);
                t -= (exp - 1) * 0.69314718055994530942;
                t *= (1 + 1e-13);
                mag_set_d(z, t);
            }
        }
        else if (fmpz_sgn(MAG_EXPREF(x)) > 0)
        {
            mag_zero(z);
        }
        else
        {
            mag_inv(z, x);
            mag_log(z, z);
        }
    }
}

void
mag_neg_log_lower(mag_t z, const mag_t x)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            mag_inf(z);
        else
            mag_zero(z);
    }
    else
    {
        double t;
        fmpz exp = MAG_EXP(x);

        if (!COEFF_IS_MPZ(exp))
        {
            if (exp >= 1)
            {
                mag_zero(z);
            }
            else if (exp > -(1000 - MAG_BITS))
            {
                t = ldexp(MAG_MAN(x), exp - MAG_BITS);
                t = -mag_d_log_upper_bound(t);
                mag_set_d_lower(z, t);
            }
            else
            {
                /* -log(2^(exp-1) * (2*v)) = -exp*log(2) - log(2*v), 1 <= 2v < 2 */
                t = MAG_MAN(x) * ldexp(1.0, 1 - MAG_BITS);
                t = -mag_d_log_upper_bound(t);
                t -= (exp - 1) * 0.69314718055994530942;
                t *= (1 - 1e-13);
                mag_set_d_lower(z, t);
            }
        }
        else if (fmpz_sgn(MAG_EXPREF(x)) > 0)
        {
            mag_zero(z);
        }
        else
        {
            mag_inv_lower(z, x);
            mag_log_lower(z, z);
        }
    }
}

void
mag_log(mag_t z, const mag_t x)
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
            if (exp <= 0 || (exp == 1 && MAG_MAN(x) == MAG_ONE_HALF))
            {
                mag_zero(z);
            }
            else if (exp < 1000)  /* safely within the double range */
            {
                t = ldexp(MAG_MAN(x), exp - MAG_BITS);
                t = mag_d_log_upper_bound(t);
                mag_set_d(z, t);
            }
            else
            {
                /* log(2^(exp-1)*(2v)) = (exp-1)*log(2) + log(2v), 1 <= 2v < 2 */
                t = MAG_MAN(x) * ldexp(1.0, 1 - MAG_BITS); /* exact */
                t = mag_d_log_upper_bound(t);
                t += (exp - 1.0) * 0.69314718055994530942;
                t *= (1 + 1e-13);  /* last step could have rounded down */
                mag_set_d(z, t);
            }
        }
        else if (fmpz_sgn(MAG_EXPREF(x)) < 0)
        {
            mag_zero(z);
        }
        else
        {
            /* log(2^exp) = exp*log(2), log(2) < 744261118/2^30 */
            fmpz_t t;
            fmpz_init(t);
            fmpz_mul_ui(t, MAG_EXPREF(x), 744261118);
            mag_set_fmpz(z, t);
            mag_mul_2exp_si(z, z, -30);
            fmpz_clear(t);
        }
    }
}

void
mag_log_lower(mag_t z, const mag_t x)
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
            if (exp <= 0 || (exp == 1 && MAG_MAN(x) == MAG_ONE_HALF))
            {
                mag_zero(z);
            }
            else if (exp < 1000)
            {
                t = ldexp(MAG_MAN(x), exp - MAG_BITS);
                t = mag_d_log_lower_bound(t);
                mag_set_d_lower(z, t);
            }
            else
            {
                /* log(2^(exp-1) * (2*v)) = exp*log(2) + log(2*v), 1 <= 2v < 2 */
                t = MAG_MAN(x) * ldexp(1.0, 1 - MAG_BITS);
                t = mag_d_log_lower_bound(t);
                t += (exp - 1) * 0.69314718055994530942;
                t *= (1 - 1e-13);
                mag_set_d_lower(z, t);
            }
        }
        else if (fmpz_sgn(MAG_EXPREF(x)) < 0)
        {
            mag_zero(z);
        }
        else
        {
            /* log(2^exp) = exp*log(2), log(2) > 744261117/2^30 */
            fmpz_t t;
            fmpz_init(t);
            fmpz_sub_ui(t, MAG_EXPREF(x), 1);  /* lower bound for x */
            fmpz_mul_ui(t, t, 744261117);
            mag_set_fmpz_lower(z, t);
            mag_mul_2exp_si(z, z, -30);
            fmpz_clear(t);
        }
    }
}

