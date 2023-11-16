/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_majorant(arb_ptr res, arb_srcptr vec, slong len, slong prec)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        arb_get_abs_ubound_arf(arb_midref(res + i), vec + i, prec);
        mag_zero(arb_radref(res + i));
    }
}

void
arb_poly_majorant(arb_poly_t res, const arb_poly_t poly, slong prec)
{
    arb_poly_fit_length(res, poly->length);
    _arb_poly_majorant(res->coeffs, poly->coeffs, poly->length, prec);
    _arb_poly_set_length(res, poly->length);
}

