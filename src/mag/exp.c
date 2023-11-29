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

static const double inverse_factorials[] = {
    1.0,
    1.0,
    0.5,
    0.16666666666666666667,
    0.041666666666666666667,
    0.0083333333333333333333,
    0.0013888888888888888889,
    0.0001984126984126984127,
    0.000024801587301587301587,
    2.7557319223985890653e-6,
    2.7557319223985890653e-7,
    2.5052108385441718775e-8,
    2.0876756987868098979e-9,
    1.6059043836821614599e-10,
    1.1470745597729724714e-11,
    7.6471637318198164759e-13
};

/* warning: requires |x| < 2^24 or thereabout */
void
_mag_exp_d(mag_t res, double x, int roundup)
{
    double u, nlog2, eps, eps2;
    slong n;

    if (roundup)
    {
        eps = 1e-13;
        eps2 = 6e-13;
    }
    else
    {
        eps = -1e-13;
        eps2 = -6e-13;
    }

    n = floor(x * 1.4426950408889634074 + 0.5);

    /* does not need to be exact */
    /* perturbed in opposite direction before subtraction */
    if (n >= 0)
        nlog2 = n * 0.69314718055994530942 * (1.0 - eps);
    else
        nlog2 = n * 0.69314718055994530942 * (1.0 + eps);

    /* perturbed (assuming |u| < 1 for absolute error to be valid) */
    u = x - nlog2 + eps;

    if (u >= -0.375 && u <= 0.375)
        u = d_polyval(inverse_factorials, 11, u) + eps2;
    else
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (roundup)
        mag_set_d(res, u);
    else
        mag_set_d_lower(res, u);

    MAG_EXP(res) += n;  /* assumes x not too large */
}

void
mag_exp_huge(mag_t res, const mag_t x)
{
    if (mag_cmp_2exp_si(x, 128) <= 0)
    {
        fmpz_t t;
        fmpz_init(t);
        mag_get_fmpz(t, x);
        MAG_MAN(res) = 729683223; /* upper bound for e */
        fmpz_set_ui(MAG_EXPREF(res), 2);
        mag_pow_fmpz(res, res, t);
        fmpz_clear(t);
    }
    else
    {
        mag_inf(res);
    }
}

void
mag_exp_huge_lower(mag_t res, const mag_t x)
{
    fmpz_t t;
    fmpz_init(t);

    if (mag_cmp_2exp_si(x, 128) <= 0)
    {
        mag_get_fmpz_lower(t, x);
    }
    else
    {
        fmpz_one(t);
        fmpz_mul_2exp(t, t, 128);
    }

    MAG_MAN(res) = 729683222; /* lower bound for e */
    fmpz_set_ui(MAG_EXPREF(res), 2);
    mag_pow_fmpz_lower(res, res, t);
    fmpz_clear(t);
}

void
mag_exp(mag_t y, const mag_t x)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            mag_one(y);
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
            MAG_MAN(y) = MAG_ONE_HALF + 1;
            fmpz_one(MAG_EXPREF(y));
        }
    }
    else
    {
        slong e = MAG_EXP(x);

        if (e <= -MAG_BITS) /* assumes MAG_BITS == 30 */
        {
            MAG_MAN(y) = MAG_ONE_HALF + 1;
            fmpz_one(MAG_EXPREF(y));
        }
        else if (e <= -(MAG_BITS / 2))  /* assumes MAG_BITS == 30 */
        {
            MAG_MAN(y) = MAG_ONE_HALF + (MAG_MAN(x) >> (1-e)) + 2;
            fmpz_one(MAG_EXPREF(y));
        }
        else if (e < 24)
        {
            _mag_exp_d(y, ldexp(MAG_MAN(x), e - MAG_BITS), 1);
        }
        else
        {
            mag_exp_huge(y, x);
        }
    }
}

void
mag_exp_lower(mag_t y, const mag_t x)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            mag_one(y);
        else
            mag_inf(y);
    }
    else if (COEFF_IS_MPZ(MAG_EXP(x)))
    {
        if (fmpz_sgn(MAG_EXPREF(x)) > 0)
            mag_exp_huge_lower(y, x);
        else
            mag_one(y);
    }
    else
    {
        slong e = MAG_EXP(x);

        if (e <= -MAG_BITS)
        {
            mag_one(y);
        }
        else if (e <= -(MAG_BITS / 2))
        {
            MAG_MAN(y) = MAG_ONE_HALF + (MAG_MAN(x) >> (1 - e));
            fmpz_one(MAG_EXPREF(y));
        }
        else if (e < 24)
        {
            _mag_exp_d(y, ldexp(MAG_MAN(x), e - MAG_BITS), 0);
        }
        else
        {
            mag_exp_huge_lower(y, x);
        }
    }
}


