/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq_mat.h"

int
fmpq_mat_solve_fmpz_mat(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
{
    if (fmpz_mat_nrows(A) <= 15)
        return fmpq_mat_solve_fmpz_mat_fraction_free(X, A, B);
    else if (2 * fmpz_mat_ncols(B) >= fmpz_mat_nrows(A))
        return fmpq_mat_solve_fmpz_mat_multi_mod(X, A, B);
    else
        return fmpq_mat_solve_fmpz_mat_dixon(X, A, B);
}

int
fmpq_mat_solve(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
{
    fmpz_mat_t Anum, Bnum;
    int success;

    if (fmpq_mat_nrows(A) <= 15)
        return fmpq_mat_solve_fraction_free(X, A, B);

    /* clear denominators rowwise and use the fmpz_mat dispatch */
    fmpz_mat_init(Anum, A->r, A->c);
    fmpz_mat_init(Bnum, B->r, B->c);
    fmpq_mat_get_fmpz_mat_rowwise_2(Anum, Bnum, NULL, A, B);
    success = fmpq_mat_solve_fmpz_mat(X, Anum, Bnum);
    fmpz_mat_clear(Anum);
    fmpz_mat_clear(Bnum);
    return success;
}
