/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "padic_poly.h"

void _padic_poly_normalise(padic_poly_t poly)
{
    slong len = poly->length;

    FMPZ_VEC_NORM(poly->coeffs, len);

    poly->length = len;
}

