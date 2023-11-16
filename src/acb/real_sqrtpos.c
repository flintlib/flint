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
acb_real_sqrtpos(acb_t res, const acb_t z, int analytic, slong prec)
{
    if (arb_is_zero(acb_imagref(z)) && !analytic)
    {
        arb_sqrtpos(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
    }
    else if (arb_is_positive(acb_realref(z)) || !arb_contains_zero(acb_imagref(z)))
    {
        acb_sqrt(res, z, prec);
    }
    else
    {
        acb_indeterminate(res);
    }
}

