/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
acb_poly_set_round(acb_poly_t dest, const acb_poly_t src, slong prec)
{
    slong len = acb_poly_length(src);

    acb_poly_fit_length(dest, len);
    _acb_vec_set_round(dest->coeffs, src->coeffs, len, prec);
    _acb_poly_set_length(dest, len);
}

