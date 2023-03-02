/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void fmpz_mod_mat_randtril(fmpz_mod_mat_t mat, flint_rand_t state, int unit)
{
    fmpz* e;
    slong i, j;

    for (i = 0; i < fmpz_mod_mat_nrows(mat); i++)
    {
        for (j = 0; j < fmpz_mod_mat_ncols(mat); j++)
        {
            e = fmpz_mod_mat_entry(mat, i, j);
            if (j < i)
            {
                fmpz_randm(e, state, mat->mod);
            }
            else if (i == j)
            {
                fmpz_randm(e, state, mat->mod);
                if (unit || fmpz_is_zero(e))
                    fmpz_one(e);
            }
            else
            {
                fmpz_zero(e);
            }
        }
    }
}

