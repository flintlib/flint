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
acb_real_max(acb_t res, const acb_t x, const acb_t y, int analytic, slong prec)
{
    arb_t t;

    if (!acb_is_finite(x) || !acb_is_finite(y))
    {
        acb_indeterminate(res);
        return;
    }

    arb_init(t);
    arb_sub(t, acb_realref(x), acb_realref(y), prec);

    if (arb_is_positive(t))
        acb_set_round(res, x, prec);
    else if (arb_is_negative(t))
        acb_set_round(res, y, prec);
    else if (!analytic)
        acb_union(res, x, y, prec);
    else
        acb_indeterminate(res);

    arb_clear(t);
}

