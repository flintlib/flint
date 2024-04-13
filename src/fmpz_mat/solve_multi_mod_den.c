/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"

int
fmpz_mat_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    int success;
    fmpq_mat_t Q;

    fmpq_mat_init(Q, fmpz_mat_nrows(X), fmpz_mat_ncols(X));
    success = fmpq_mat_solve_fmpz_mat_multi_mod(Q, A, B);

    if (success)
        fmpq_mat_get_fmpz_mat_matwise(X, den, Q);

    fmpq_mat_clear(Q);
    return success;
}
