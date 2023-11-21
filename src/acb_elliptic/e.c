/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"

void
acb_elliptic_e(acb_t res, const acb_t m, slong prec)
{
    if (acb_is_zero(m))
    {
        acb_const_pi(res, prec);
        acb_mul_2exp_si(res, res, -1);
    }
    else if (acb_is_one(m))
    {
        acb_one(res);
    }
    else
    {
        acb_struct t[2];

        acb_init(t + 0);
        acb_init(t + 1);

        acb_elliptic_k_jet(t, m, 2, prec);
        acb_mul(t + 1, t + 1, m, prec);
        acb_mul_2exp_si(t + 1, t + 1, 1);
        acb_add(t, t, t + 1, prec);
        acb_sub_ui(t + 1, m, 1, prec);
        acb_mul(res, t, t + 1, prec);
        acb_neg(res, res);

        acb_clear(t + 0);
        acb_clear(t + 1);
    }
}

