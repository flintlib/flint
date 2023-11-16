/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_acos(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_one(z))
    {
        acb_zero(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_one(t);
        if (arb_is_zero(acb_imagref(z))
                && arb_gt(acb_realref(z), acb_realref(t)))
        {
            /* pure imaginary on (1,inf) */
            acb_asin(res, z, prec);
            acb_neg(res, res);
            arb_zero(acb_realref(res));
        }
        else
        {
            acb_asin(res, z, prec);
            acb_const_pi(t, prec);
            acb_mul_2exp_si(t, t, -1);
            acb_sub(res, t, res, prec);
        }
        acb_clear(t);
    }
}

