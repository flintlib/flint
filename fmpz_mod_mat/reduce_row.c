/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"

slong _fmpz_mod_mat_reduce_row(fmpz_mod_mat_t A, slong * P, slong * L, 
                                             slong m, const fmpz_mod_ctx_t ctx)
{
    slong n = fmpz_mod_mat_ncols(A), i, j, r;
    fmpz_t h;

    fmpz_init(h);

    for (i = 0; i < n; i++)
    {
        if (!fmpz_is_zero(fmpz_mod_mat_entry(A, m, i)))
        {
            r = P[i];
            if (r != -WORD(1))
            {
                for (j = i + 1; j < L[r]; j++)
                {
                    fmpz_mod_mul(h, fmpz_mod_mat_entry(A, r, j),
                                    fmpz_mod_mat_entry(A, m, i), ctx);
                    fmpz_mod_sub(fmpz_mod_mat_entry(A, m, j),
                                 fmpz_mod_mat_entry(A, m, j), h, ctx);
                }
 
                fmpz_zero(fmpz_mod_mat_entry(A, m, i));
            }
            else
            {
                fmpz_mod_inv(h, fmpz_mod_mat_entry(A, m, i), ctx);
                fmpz_one(fmpz_mod_mat_entry(A, m, i));
           
                for (j = i + 1; j < L[m]; j++)
                    fmpz_mod_mul(fmpz_mod_mat_entry(A, m, j),
                                 fmpz_mod_mat_entry(A, m, j), h, ctx);
               
                P[i] = m;

                fmpz_clear(h);
                return i;
            }
        }
    }

    fmpz_clear(h);
    return -WORD(1);
}

slong fmpz_mat_reduce_row(fmpz_mod_mat_t A, slong * P, slong * L, slong m)
{
    slong res;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_ctx_init(ctx, A->mod);
    res = _fmpz_mod_mat_reduce_row(A, P, L, m, ctx);
    fmpz_mod_ctx_clear(ctx);
    return res;
}

