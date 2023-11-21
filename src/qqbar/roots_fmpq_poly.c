/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_roots_fmpq_poly(qqbar_ptr res, const fmpq_poly_t poly, int flags)
{
    fmpz_poly_t t; /* fake an fmpz_poly */
    t->coeffs = poly->coeffs;
    t->length = poly->length;
    t->alloc = poly->alloc;
    qqbar_roots_fmpz_poly(res, t, flags);
}
