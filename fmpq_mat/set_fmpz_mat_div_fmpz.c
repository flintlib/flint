/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void
fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t X, const fmpz_mat_t Xnum,
                                                const fmpz_t den)
{
    slong i, j;

    if (fmpz_is_one(den))
    {
        fmpq_mat_set_fmpz_mat(X, Xnum);
    }
    else if (*den == WORD(-1))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set(t, den);

        for (i = 0; i < Xnum->r; i++)
        {
            for (j = 0; j < Xnum->c; j++)
            {
                fmpz_neg(fmpq_mat_entry_num(X, i, j),
                    fmpz_mat_entry(Xnum, i, j));
                fmpz_one(fmpq_mat_entry_den(X, i, j));
            }
        }

        fmpz_clear(t);
    }
    else
    {
        for (i = 0; i < Xnum->r; i++)
        {
            for (j = 0; j < Xnum->c; j++)
            {
                fmpz_set(fmpq_mat_entry_num(X, i, j),
                                fmpz_mat_entry(Xnum, i, j));
                fmpz_set(fmpq_mat_entry_den(X, i, j), den);
                fmpq_canonicalise(fmpq_mat_entry(X, i, j));
            }
        }
    }
}
