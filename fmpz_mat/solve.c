/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "perm.h"

int
fmpz_mat_solve(fmpz_mat_t X, fmpz_t den,
                    const fmpz_mat_t A, const fmpz_mat_t B)
{
    if (fmpz_mat_nrows(A) <= 3)
        return fmpz_mat_solve_cramer(X, den, A, B);
    else if (fmpz_mat_nrows(A) <= 15)
        return fmpz_mat_solve_fflu(X, den, A, B);
    else if (fmpz_mat_ncols(B) == 1)
        return fmpz_mat_solve_dixon_den(X, den, A, B);
    else
        return fmpz_mat_solve_multi_mod_den(X, den, A, B);
}
