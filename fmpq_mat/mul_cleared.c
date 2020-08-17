/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void fmpq_mat_mul_cleared(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
{
    slong i, j;

    fmpz_mat_t Aclear;
    fmpz_mat_t Bclear;
    fmpz_mat_t Cclear;

    fmpz * Aden;
    fmpz * Bden;

    fmpz_mat_init(Aclear, A->r, A->c);
    fmpz_mat_init(Bclear, B->r, B->c);
    fmpz_mat_init(Cclear, A->r, B->c);

    Aden = _fmpz_vec_init(A->r);
    Bden = _fmpz_vec_init(B->c);

    fmpq_mat_get_fmpz_mat_rowwise(Aclear, Aden, A);
    fmpq_mat_get_fmpz_mat_colwise(Bclear, Bden, B);

    fmpz_mat_mul(Cclear, Aclear, Bclear);

    for (i = 0; i < C->r; i++)
    {
        for (j = 0; j < C->c; j++)
        {
            fmpz_set(fmpq_mat_entry_num(C, i, j), fmpz_mat_entry(Cclear, i, j));
            fmpz_mul(fmpq_mat_entry_den(C, i, j), Aden + i, Bden + j);
            fmpq_canonicalise(fmpq_mat_entry(C, i, j));
        }
    }

    fmpz_mat_clear(Aclear);
    fmpz_mat_clear(Bclear);
    fmpz_mat_clear(Cclear);

    _fmpz_vec_clear(Aden, A->r);
    _fmpz_vec_clear(Bden, B->c);
}
