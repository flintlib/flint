/*
    Copyright (C) 2005-2009 Damien Stehle
    Copyright (C) 2007 David Cade
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_randajtai(fmpz_mat_t mat, flint_rand_t state, double alpha)
{
    const slong c = mat->c, r = mat->r, d = r;

    slong i, j;
    fmpz_t tmp;

    if (c != r)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_ajtai): Non-square matrix.\n");
    }

    fmpz_init(tmp);

    for (i = 0; i < d; i++)
    {
        flint_bitcnt_t bits = (flint_bitcnt_t) pow((double) (2 * d - i), alpha);

        fmpz_one(tmp);
        fmpz_mul_2exp(tmp, tmp, bits);
        fmpz_sub_ui(tmp, tmp, 1);
        fmpz_randm(fmpz_mat_entry(mat, i, i), state, tmp);
        fmpz_add_ui(fmpz_mat_entry(mat, i, i), fmpz_mat_entry(mat, i, i), 2);
        fmpz_fdiv_q_2exp(fmpz_mat_entry(mat, i, i), fmpz_mat_entry(mat, i, i), 1);

        for (j = i + 1; j < d; j++)
        {
            fmpz_randm(fmpz_mat_entry(mat, j, i), state, tmp);
            if (n_randint(state, 2))
                fmpz_neg(fmpz_mat_entry(mat, j, i), fmpz_mat_entry(mat, j, i));
            fmpz_zero(fmpz_mat_entry(mat, i, j));
        }
    }

    fmpz_clear(tmp);
}
