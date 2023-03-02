/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_mul(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
    const fmpz_poly_mat_t B)
{
    if (A->r < 8 || B->r < 8 || B->c < 8)
    {
        fmpz_poly_mat_mul_classical(C, A, B);
    }
    else
    {
        fmpz_poly_mat_mul_KS(C, A, B);
    }
}
