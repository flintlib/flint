/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_real_abs(acb_t res, const acb_t z, int analytic, slong prec)
{
    if (!acb_is_finite(z) || (analytic && arb_contains_zero(acb_realref(z))))
    {
        acb_indeterminate(res);
    }
    else
    {
        if (arb_is_nonnegative(acb_realref(z)))
        {
            acb_set_round(res, z, prec);
        }
        else if (arb_is_negative(acb_realref(z)))
        {
            acb_neg_round(res, z, prec);
        }
        else
        {
            acb_t t;
            acb_init(t);
            acb_neg(t, z);
            acb_union(res, z, t, prec);
            acb_clear(t);
        }
    }
}
