/*
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2)
{
    FLINT_SWAP(fmpz_poly_struct, *poly1, *poly2);
}
