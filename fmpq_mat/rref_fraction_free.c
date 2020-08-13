/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

slong
fmpq_mat_rref_fraction_free(fmpq_mat_t B, const fmpq_mat_t A)
{
    fmpz_mat_t Aclear;
    fmpz_t den;
    slong rank;

    if (fmpq_mat_is_empty(A))
        return 0;

    fmpz_mat_init(Aclear, A->r, A->c);
    fmpq_mat_get_fmpz_mat_rowwise(Aclear, NULL, A);
    fmpz_init(den);

    rank = fmpz_mat_rref(Aclear, den, Aclear);

    if (rank == 0)
        fmpq_mat_zero(B);
    else
        fmpq_mat_set_fmpz_mat_div_fmpz(B, Aclear, den);

    fmpz_mat_clear(Aclear);
    fmpz_clear(den);

    return rank;
}
