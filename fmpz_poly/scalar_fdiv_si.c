/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
fmpz_poly_scalar_fdiv_si(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                         slong x)
{
    if (x == 0)
    {
        flint_printf("Exception (fmpz_poly_scalar_fdiv_si). Division by zero.\n");
        flint_abort();
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_fdiv_q_si(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}
