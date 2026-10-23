/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

int
fmpq_mat_can_solve(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
{
    if (fmpq_mat_nrows(A) <= 15)
        return fmpq_mat_can_solve_fraction_free(X, A, B);
    else
    {
        fmpz_mat_t Anum, Bnum, Z;
        fmpz_t den;
        int success;
        fmpz_mat_init(Anum, A->r, A->c);
        fmpz_mat_init(Bnum, B->r, B->c);
        fmpz_mat_init(Z, A->c, B->c);
        fmpz_init(den);
        fmpq_mat_get_fmpz_mat_rowwise_2(Anum, Bnum, NULL, A, B);
        success = _fmpz_mat_can_solve_auto_den(Z, den, Anum, Bnum);
        if (success)
            fmpq_mat_set_fmpz_mat_div_fmpz(X, Z, den);
        fmpz_mat_clear(Anum);
        fmpz_mat_clear(Bnum);
        fmpz_mat_clear(Z);
        fmpz_clear(den);
        return success;
    }
}
