/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "ulong_extras.h"

void
fmpz_poly_mat_randtest_sparse(fmpz_poly_mat_t A, flint_rand_t state, slong len,
    flint_bitcnt_t bits, float density)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (n_randint(state, 1000) < density * 1000)
            {
                slong l = n_randint(state, len + 1);
                l = FLINT_MAX(l, 1);
                fmpz_poly_randtest_not_zero(fmpz_poly_mat_entry(A, i, j),
                    state, l, bits);
            }
            else
            {
                fmpz_poly_zero(fmpz_poly_mat_entry(A, i, j));
            }
        }
    }
}
