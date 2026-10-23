/*
    Copyright (C) 2020 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"

int
fmpq_mat_can_solve_fmpz_mat_multi_mod(fmpq_mat_t X,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    fmpz_mat_t Z;
    fmpz_t den;
    int success;

    if (A->r != B->r || A->c != X->r || X->c != B->c)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_mat_can_solve_fmpz_mat_multi_mod). Incompatible matrix dimensions.\n");
    }

    fmpz_mat_init(Z, A->c, B->c);
    fmpz_init(den);

    success = fmpz_mat_can_solve_multi_mod_den(Z, den, A, B);
    if (success)
        fmpq_mat_set_fmpz_mat_div_fmpz(X, Z, den);

    fmpz_mat_clear(Z);
    fmpz_clear(den);
    return success;
}

int
fmpq_mat_can_solve_multi_mod(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
{
    fmpz_mat_t Anum;
    fmpz_mat_t Bnum;
    int success;

    if (A->r != B->r || A->c != X->r || X->c != B->c)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_mat_can_solve_multi_mod). Incompatible matrix dimensions.\n");
    }

    if (A->r == 0 || B->c == 0)
    {
        fmpq_mat_zero(X);
        return 1;
    }

    if (A->c == 0)
    {
        fmpq_mat_zero(X);
        return fmpq_mat_is_zero(B);
    }

    fmpz_mat_init(Anum, A->r, A->c);
    fmpz_mat_init(Bnum, B->r, B->c);

    fmpq_mat_get_fmpz_mat_rowwise_2(Anum, Bnum, NULL, A, B);
    success = fmpq_mat_can_solve_fmpz_mat_multi_mod(X, Anum, Bnum);

    fmpz_mat_clear(Anum);
    fmpz_mat_clear(Bnum);

    return success;
}
