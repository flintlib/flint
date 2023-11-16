/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_modular.h"
#include "acb_elliptic.h"

void
_acb_poly_elliptic_k_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
{
    _acb_elliptic_k_series(res, z, zlen, len, prec);
}

void
acb_poly_elliptic_k_series(acb_poly_t res, const acb_poly_t z, slong n, slong prec)
{
    acb_elliptic_k_series(res, z, n, prec);
}
