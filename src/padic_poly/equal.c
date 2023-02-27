/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

int padic_poly_equal(const padic_poly_t f, const padic_poly_t g)
{
    if (f == g)
    {
        return 1;
    }

    if (f->length != g->length || f->val != g->val)
    {
        return 0;
    }

    return _fmpz_vec_equal(f->coeffs, g->coeffs, f->length);
}

