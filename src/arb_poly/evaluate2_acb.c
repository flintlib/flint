/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_evaluate2_acb(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec)
{
    _arb_poly_evaluate2_acb_rectangular(y, z, f, len, x, prec);
}

void
arb_poly_evaluate2_acb(acb_t r, acb_t s, const arb_poly_t f, const acb_t a, slong prec)
{
    _arb_poly_evaluate2_acb(r, s, f->coeffs, f->length, a, prec);
}

