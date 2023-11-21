/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_sqrts(acb_t y1, acb_t y2, const acb_t x, slong prec)
{
    if (acb_contains_zero(x))
    {
        acb_sqrt(y1, x, prec);
        acb_neg(y2, y1);
        acb_union(y1, y1, y2, prec);
        acb_set(y2, y1);
    }
    else if (arb_contains_zero(acb_imagref(x)) && arb_is_negative(acb_realref(x)))
    {
        acb_neg(y1, x);
        acb_sqrt(y1, y1, prec);
        acb_mul_onei(y1, y1);
        acb_neg(y2, y1);
    }
    else
    {
        acb_sqrt(y1, x, prec);
        acb_neg(y2, y1);
    }
}
