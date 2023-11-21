/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_expint(acb_t res, const acb_t s, const acb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_sub_ui(t, s, 1, prec);
    acb_neg(t, t);
    acb_hypgeom_gamma_upper(res, t, z, 2, prec);
    acb_clear(t);
}

