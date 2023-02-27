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
fmpz_mat_CRT_ui(fmpz_mat_t res, const fmpz_mat_t mat1,
                        const fmpz_t m1, const nmod_mat_t mat2, int sign)
{
    slong i, j;
    mp_limb_t c;
    mp_limb_t m2 = mat2->mod.n;
    mp_limb_t m2inv = mat2->mod.ninv;
    fmpz_t m1m2;

    c = fmpz_fdiv_ui(m1, m2);
    c = n_invmod(c, m2);

    if (c == 0)
    {
        flint_printf("Exception (fmpz_mat_CRT_ui). m1 not invertible modulo m2.\n");
        flint_abort();
    }

    fmpz_init(m1m2);
    fmpz_mul_ui(m1m2, m1, m2);

    for (i = 0; i < mat1->r; i++)
    {
        for (j = 0; j < mat1->c; j++)
            _fmpz_CRT_ui_precomp(fmpz_mat_entry(res, i, j),
                    fmpz_mat_entry(mat1, i, j), m1,
                    nmod_mat_entry(mat2, i, j), m2, m2inv, m1m2, c, sign);
    }

    fmpz_clear(m1m2);
}

