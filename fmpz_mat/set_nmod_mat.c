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
fmpz_mat_set_nmod_mat(fmpz_mat_t A, const nmod_mat_t Amod)
{
    slong i, j;

    for (i = 0; i < Amod->r; i++)
        for (j = 0; j < Amod->c; j++)
            fmpz_set_ui_smod(fmpz_mat_entry(A, i, j),
                             nmod_mat_entry(Amod, i, j), Amod->mod.n);
}
