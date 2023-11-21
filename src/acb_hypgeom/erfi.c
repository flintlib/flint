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
acb_hypgeom_erfi(acb_t res, const acb_t z, slong prec)
{
    acb_mul_onei(res, z);
    acb_hypgeom_erf(res, res, prec);
    acb_mul_onei(res, res);
    acb_neg(res, res);
}

