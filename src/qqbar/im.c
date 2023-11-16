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
qqbar_im(qqbar_t res, const qqbar_t x)
{
    if (qqbar_sgn_im(x) == 0)
    {
        qqbar_zero(res);
    }
    else
    {
        qqbar_t t;
        qqbar_init(t);

        if (qqbar_sgn_re(x) == 0)
        {
            qqbar_i(t);
            qqbar_mul(res, x, t);
            qqbar_neg(res, res);
        }
        else
        {
            qqbar_conj(t, x);
            qqbar_sub(res, x, t);
            qqbar_i(t);
            qqbar_mul(res, res, t);
            qqbar_neg(res, res);
            qqbar_mul_2exp_si(res, res, -1);
        }

        arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));
        qqbar_clear(t);
    }
}
