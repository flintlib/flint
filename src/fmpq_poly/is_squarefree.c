/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpq_poly.h"

int fmpq_poly_is_squarefree(const fmpq_poly_t poly)
{
    return _fmpz_poly_is_squarefree(poly->coeffs, poly->length);
}
