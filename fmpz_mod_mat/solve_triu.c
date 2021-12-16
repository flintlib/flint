/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void fmpz_mod_mat_solve_triu(fmpz_mod_mat_t X, const fmpz_mod_mat_t U,
                             const fmpz_mod_mat_t B, int unit)
{
    if (B->mat->r < FMPZ_MOD_MAT_SOLVE_TRI_ROWS_CUTOFF ||
        B->mat->c < FMPZ_MOD_MAT_SOLVE_TRI_COLS_CUTOFF)
    {
        fmpz_mod_mat_solve_triu_classical(X, U, B, unit);
    }
    else
    {
        fmpz_mod_mat_solve_triu_recursive(X, U, B, unit);
    }
}

