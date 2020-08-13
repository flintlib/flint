/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_scalar_addmul_nmod_mat_fmpz(fmpz_mat_t B,
                        const nmod_mat_t A, const fmpz_t c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_addmul_ui(fmpz_mat_entry(B,i,j), c, nmod_mat_entry(A,i,j));
}
