/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
fmpz_mat_set(fmpz_mat_t mat1, const fmpz_mat_t mat2)
{
    if (mat1 != mat2)
    {
        slong i;

        if (mat2->r && mat2->c)
            for (i = 0; i < mat2->r; i++)
                _fmpz_vec_set(fmpz_mat_row(mat1, i), fmpz_mat_row(mat2, i), mat2->c);
    }
}

void
fmpz_mat_set_nmod_mat(fmpz_mat_t A, const nmod_mat_t Amod)
{
    slong i, j;

    for (i = 0; i < Amod->r; i++)
        for (j = 0; j < Amod->c; j++)
            fmpz_set_ui_smod(fmpz_mat_entry(A, i, j),
                             nmod_mat_entry(Amod, i, j), Amod->mod.n);
}

void
fmpz_mat_set_nmod_mat_unsigned(fmpz_mat_t A, const nmod_mat_t Amod)
{
    slong i, j;

    for (i = 0; i < Amod->r; i++)
        for (j = 0; j < Amod->c; j++)
            fmpz_set_ui(fmpz_mat_entry(A, i, j), nmod_mat_entry(Amod, i, j));
}
