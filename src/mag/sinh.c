/*
    Copyright (C) 2014, 2018 Fredrik Johansson

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
    0.16666666666666666667,
    0.0083333333333333333333,
    0.0001984126984126984127,
    2.7557319223985890653e-6,
    2.5052108385441718775e-8,
    1.6059043836821614599e-10,
    7.6471637318198164759e-13
};

static double _sinh(double x)
{
    x = x * (1.0 / 27.0);
    x = x * d_polyval(inverse_factorials, 8, x * x);
    x = x * (3.0 + 4.0 * (x * x));
    x = x * (3.0 + 4.0 * (x * x));
    x = x * (3.0 + 4.0 * (x * x));
    return x;
}

void
mag_sinh(mag_t res, const mag_t x)
{
    if (mag_is_special(x))
    {
        mag_set(res, x);
    }
    else
    {
        if (mag_cmp_2exp_si(x, -30) < 0)
        {
            mag_expm1(res, x);
        }
        else if (mag_cmp_2exp_si(x, 4) > 0)
        {
            mag_exp(res, x);
            mag_mul_2exp_si(res, res, -1);
        }
        else
        {
            double v;
            v = mag_get_d(x);
            v = _sinh(v) * (1 + 1e-12);
            mag_set_d(res, v);
        }
    }
}

void
mag_sinh_lower(mag_t res, const mag_t x)
{
    if (mag_is_special(x))
    {
        mag_set(res, x);
    }
    else
    {
        if (mag_cmp_2exp_si(x, -15) < 0)
        {
            mag_set(res, x);
        }
        else if (mag_cmp_2exp_si(x, 4) > 0)
        {
            mag_t t;
            mag_init(t);
            mag_exp_lower(t, x);
            mag_expinv(res, x);
            mag_sub(res, t, res);
            mag_mul_2exp_si(res, res, -1);
            mag_clear(t);
        }
        else
        {
            double v;
            v = mag_get_d(x);
            v = _sinh(v) * (1 - 1e-12);
            mag_set_d_lower(res, v);
        }
    }
}
