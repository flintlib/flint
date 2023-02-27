/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "perm.h"

int
fmpz_poly_mat_solve(fmpz_poly_mat_t X, fmpz_poly_t den,
                    const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
{
    return fmpz_poly_mat_solve_fflu(X, den, A, B);
}
