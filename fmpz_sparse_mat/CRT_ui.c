/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_sparse_mat.h"
#include "fmpz_sparse_mat.h"

void
fmpz_sparse_mat_CRT_ui(fmpz_sparse_mat_t res, const fmpz_sparse_mat_t mat1,
                        const fmpz_t m1, const nmod_sparse_mat_t mat2, int sign)
{
    slong i;
    mp_limb_t m1i_m2 = n_invmod(fmpz_fdiv_ui(m1, mat2->mod.n), mat2->mod.n);

    if (m1i_m2 == 0)
    {
        flint_printf("Exception (fmpz_mat_CRT_ui). m1 not invertible modulo m2.\n");
        flint_abort();
    }
    for (i = 0; i < mat1->r; i++)
        fmpz_sparse_vec_CRT_ui(&res->rows[i], &mat1->rows[i], m1, &mat2->rows[i], mat2->mod, m1i_m2, sign);
}

