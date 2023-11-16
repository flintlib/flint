/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
arb_poly_set_round(arb_poly_t dest, const arb_poly_t src, slong prec)
{
    slong len = arb_poly_length(src);

    arb_poly_fit_length(dest, len);
    _arb_vec_set_round(dest->coeffs, src->coeffs, len, prec);
    _arb_poly_set_length(dest, len);
}

