/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

int fmpz_mod_mat_reduce_row(slong * column, fmpz_mod_mat_t A, slong * P, slong * L,
                                         slong m, const fmpz_mod_ctx_t ctx)
{
    slong n = fmpz_mod_mat_ncols(A, ctx), i, j, r, res = -WORD(1);
    fmpz_t d, h;
    int status = GR_SUCCESS;

    fmpz_init(d);
    fmpz_init(h);

    for (i = 0; i < n; i++)
    {
        if (i != 0)
            fmpz_mod_set_fmpz(fmpz_mod_mat_entry(A, m, i),
                              fmpz_mod_mat_entry(A, m, i), ctx);

        if (!fmpz_is_zero(fmpz_mod_mat_entry(A, m, i)))
        {
            r = P[i];

            if (r != -WORD(1))
            {
                for (j = i + 1; j < L[r]; j++)
                {
                    fmpz_submul(fmpz_mod_mat_entry(A, m, j),
                        fmpz_mod_mat_entry(A, r, j),
                        fmpz_mod_mat_entry(A, m, i));
                }

                fmpz_zero(fmpz_mod_mat_entry(A, m, i));
            }
            else
            {
                /* fmpz_mod_inv(h, fmpz_mod_mat_entry(A, m, i), ctx); */
                fmpz_gcdinv(d, h, fmpz_mod_mat_entry(A, m, i), ctx->n);
                if (!fmpz_is_one(d))
                {
                    status = GR_DOMAIN;
                    goto cleanup;
                }

                fmpz_mod_set_ui(fmpz_mod_mat_entry(A, m, i), 1, ctx);

                for (j = i + 1; j < L[m]; j++)
                {
                    fmpz_mod_set_fmpz(fmpz_mod_mat_entry(A, m, j),
                        fmpz_mod_mat_entry(A, m, j), ctx);
                    fmpz_mod_mul(fmpz_mod_mat_entry(A, m, j),
                        fmpz_mod_mat_entry(A, m, j), h, ctx);
                }

                P[i] = m;
                res = i;
                break;

            }
        }
    }

cleanup:
    fmpz_clear(h);
    fmpz_clear(d);
    *column = res;
    return status;
}
