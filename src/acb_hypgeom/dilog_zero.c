/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_dilog_zero(acb_t res, const acb_t z, slong prec)
{
    if (prec < 40000 ||
        (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -prec / 1000) < 0
      && arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -prec / 1000) < 0))
    {
        acb_hypgeom_dilog_zero_taylor(res, z, prec);
    }
    else
    {
        acb_t z0;
        acb_init(z0);
        acb_hypgeom_dilog_bitburst(res, z0, z, prec);
        acb_hypgeom_dilog_zero_taylor(z0, z0, prec);
        acb_add(res, res, z0, prec);
        acb_clear(z0);
    }
}

