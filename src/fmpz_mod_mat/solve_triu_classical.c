/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"

void fmpz_mod_mat_solve_triu_classical(fmpz_mod_mat_t X,
                                       const fmpz_mod_mat_t U,
                                       const fmpz_mod_mat_t B,
                                       int unit)
{
    slong i, j, n, m;
    fmpz * inv, * tmp;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init(ctx, U->mod);

    n = fmpz_mod_mat_nrows(U);
    m = fmpz_mod_mat_ncols(B);

    if (!unit)
    {
        inv = _fmpz_vec_init(n);
        for (i = 0; i < n; i++)
            fmpz_mod_inv(inv + i, fmpz_mod_mat_entry(U, i, i), ctx);
    }
    else
        inv = NULL;

    tmp = _fmpz_vec_init(n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            fmpz_set(tmp + j, fmpz_mod_mat_entry(X, j, i));

        for (j = n - 1; j >= 0; j--)
        {
            fmpz_t s;
            fmpz_init(s);
            _fmpz_mod_vec_dot(s, U->mat->rows[j] + j + 1, tmp + j + 1,
                                   n - j - 1, ctx);
            fmpz_mod_sub(s, fmpz_mod_mat_entry(B, j, i), s, ctx);
            if (!unit)
                fmpz_mod_mul(s, s, inv + j, ctx);
            fmpz_set(tmp + j, s);
            fmpz_clear(s);
        }

        for (j = 0; j < n; j++)
            fmpz_mod_mat_set_entry(X, j, i, tmp + j);
    }

    _fmpz_vec_clear(tmp, n);
    if (!unit)
        _fmpz_vec_clear(inv, n);

    fmpz_mod_ctx_clear(ctx);
}

