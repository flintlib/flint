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
qqbar_re(qqbar_t res, const qqbar_t x)
{
    if (qqbar_sgn_im(x) == 0)
    {
        qqbar_set(res, x);
    }
    else if (qqbar_sgn_re(x) == 0)
    {
        qqbar_zero(res);
    }
    else
    {
        qqbar_t t;
        qqbar_init(t);
        qqbar_conj(t, x);
        qqbar_add(res, x, t);
        arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));
        qqbar_mul_2exp_si(res, res, -1);
        qqbar_clear(t);
    }
}
