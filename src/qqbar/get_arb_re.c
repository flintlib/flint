/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_get_arb_re(arb_t res, const qqbar_t x, slong prec)
{
    if (qqbar_sgn_re(x) == 0)
    {
        arb_zero(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        qqbar_get_acb(t, x, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}
