/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_poly_erf_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
{
    _acb_hypgeom_erf_series(g, h, hlen, n, prec);
}

void
acb_poly_erf_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
{
    acb_hypgeom_erf_series(g, h, n, prec);
}

