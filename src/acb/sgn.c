/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_sgn(acb_t res, const acb_t z, slong prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_sgn(acb_realref(res), acb_realref(z));
        arb_zero(acb_imagref(res));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_sgn(acb_imagref(res), acb_imagref(z));
        arb_zero(acb_realref(res));
    }
    else
    {
        arb_t t;
        arb_init(t);
        acb_abs(t, z, prec);
        arb_inv(t, t, prec);

        if (arb_is_finite(t))
        {
            acb_mul_arb(res, z, t, prec);
        }
        else
        {
            arf_zero(arb_midref(acb_realref(res)));
            mag_one(arb_radref(acb_realref(res)));
            arb_set(acb_imagref(res), acb_realref(res));
        }

        arb_clear(t);
    }
}

