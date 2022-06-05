/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_mat-impl.h"


void
fmpz_poly_mat_init_set(fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
{
    fmpz_poly_mat_init(A, B->r, B->c);
    fmpz_poly_mat_set(A, B);
}
