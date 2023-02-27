/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"

void
nmod_mat_solve_triu(nmod_mat_t X, const nmod_mat_t U,
                                    const nmod_mat_t B, int unit)
{
    if (B->r < NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF ||
        B->c < NMOD_MAT_SOLVE_TRI_COLS_CUTOFF)
    {
        nmod_mat_solve_triu_classical(X, U, B, unit);
    }
    else
    {
        nmod_mat_solve_triu_recursive(X, U, B, unit);
    }
}
