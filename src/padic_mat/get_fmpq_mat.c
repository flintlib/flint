/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_mat.h"

void padic_mat_get_fmpq_mat(fmpq_mat_t B, 
                            const padic_mat_t A, const padic_ctx_t ctx)
{
    if (!padic_mat_is_empty(A))
    {
        if (padic_mat_is_zero(A))
        {
            fmpq_mat_zero(B);
        }
        else
        {
            fmpz_t f;
            slong i, j;

            fmpz_init(f);
            fmpz_pow_ui(f, ctx->p, FLINT_ABS(padic_mat_val(A)));
            for (i = 0; i < B->r; i++)
                for (j = 0; j < B->c; j++)
                {
                    if (padic_mat_val(A) >= 0)
                    {
                        fmpz_mul(fmpq_mat_entry_num(B, i, j), padic_mat_entry(A, i, j), f);
                        fmpz_one(fmpq_mat_entry_den(B, i, j));
                    }
                    else
                    {
                        fmpz_set(fmpq_mat_entry_num(B, i, j), padic_mat_entry(A, i, j));
                        fmpz_set(fmpq_mat_entry_den(B, i, j), f);
                        fmpq_canonicalise(fmpq_mat_entry(B, i, j));
                    }
                }
            fmpz_clear(f);
        }
    }
}

