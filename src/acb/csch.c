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
acb_csch(acb_t res, const acb_t z, slong prec)
{
    if (acb_contains_zero(z) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
    }
    else if (arb_is_zero(acb_imagref(z)))
    {
        arb_csch(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_csc(acb_imagref(res), acb_imagref(z), prec);
        arb_neg(acb_imagref(res), acb_imagref(res));
        arb_zero(acb_realref(res));
    }
    else
    {
        if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) > 0)
        {
            acb_t t;
            acb_init(t);

            if (arf_sgn(arb_midref(acb_realref(z))) > 0)
            {
                acb_neg(t, z);
                acb_exp(t, t, prec + 4);
                acb_mul(res, t, t, prec + 4);
                acb_sub_ui(res, res, 1, prec + 4);
                acb_div(res, t, res, prec);
                acb_neg(res, res);
            }
            else
            {
                acb_exp(t, z, prec + 4);
                acb_mul(res, t, t, prec + 4);
                acb_sub_ui(res, res, 1, prec + 4);
                acb_div(res, t, res, prec);
            }

            acb_mul_2exp_si(res, res, 1);
            acb_clear(t);
        }
        else
        {
            acb_sinh(res, z, prec + 4);
            acb_inv(res, res, prec);
        }
    }
}

