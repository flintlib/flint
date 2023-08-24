/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_sqrts(acb_t y1, acb_t y2, const acb_t x, slong prec)
{
    if (arb_contains_zero(acb_imagref(x)) && arb_is_negative(acb_realref(x)))
    {
        acb_neg(y1, x);
        acb_sqrt(y1, y1, prec);
        acb_mul_onei(y1, y1);
    }
    else
    {
        acb_sqrt(y1, x, prec);
    }
    acb_neg(y2, y1);
}
