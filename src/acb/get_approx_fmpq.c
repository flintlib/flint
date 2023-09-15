/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

int acb_get_approx_fmpq(fmpq_t y, const acb_t x, slong prec)
{
    int res = 0;
    mag_t im;

    mag_init(im);
    arb_get_mag(im, acb_imagref(x));

    if (mag_cmp_2exp_si(im, -prec / 8) < 0)
    {
        res = arb_get_approx_fmpq(y, acb_realref(x), prec);
    }

    mag_clear(im);
    return res;
}
