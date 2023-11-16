/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_cosh(mag_t res, const mag_t x)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
    }
    else if (mag_is_inf(x))
    {
        mag_inf(res);
    }
    else
    {
        mag_t t;
        mag_init(t);
        mag_exp(t, x);
        mag_expinv(res, x);
        mag_add(res, res, t);
        mag_mul_2exp_si(res, res, -1);
        mag_clear(t);
    }
}

void
mag_cosh_lower(mag_t res, const mag_t x)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
    }
    else if (mag_is_inf(x))
    {
        mag_inf(res);
    }
    else
    {
        mag_t t;
        mag_init(t);
        mag_exp_lower(t, x);
        mag_expinv_lower(res, x);
        mag_add_lower(res, res, t);
        mag_mul_2exp_si(res, res, -1);
        mag_clear(t);
    }
}
