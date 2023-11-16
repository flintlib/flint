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
acb_real_heaviside(acb_t res, const acb_t z, int analytic, slong prec)
{
    acb_real_sgn(res, z, analytic, prec);

    if (acb_is_finite(res))
    {
        acb_add_ui(res, res, 1, prec);
        acb_mul_2exp_si(res, res, -1);
    }
}

