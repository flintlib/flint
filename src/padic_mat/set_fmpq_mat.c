/*
    Copyright (C) 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <limits.h>

#include "padic_mat.h"

void padic_mat_set_fmpq_mat(padic_mat_t B, 
                            const fmpq_mat_t A, const padic_ctx_t ctx)
{
    if (!fmpq_mat_is_empty(A))
    {
        const slong N = padic_mat_prec(B);

        slong i, j, m = WORD_MAX, v, w;
        fmpz_t f, g, s, t;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(s);
        fmpz_init(t);

        /* Find min valuation m */
        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                if (!fmpq_is_zero(fmpq_mat_entry(A, i, j)))
                {
                    v = fmpz_remove(t, fmpq_mat_entry_num(A, i, j), ctx->p);
                    w = fmpz_remove(t, fmpq_mat_entry_den(A, i, j), ctx->p);
                    m = FLINT_MIN(m, v - w);
                }

        if (m >= N)
        {
            padic_mat_zero(B);
        }
        else
        {
            for (i = 0; i < A->r; i++)
                for (j = 0; j < A->c; j++)
                    if (fmpq_is_zero(fmpq_mat_entry(A, i, j)))
                    {
                        fmpz_zero(padic_mat_entry(B, i, j));
                    }
                    else
                    {
                        v = fmpz_remove(s, fmpq_mat_entry_num(A, i, j), ctx->p);
                        w = fmpz_remove(t, fmpq_mat_entry_den(A, i, j), ctx->p);
                        if (v - w >= N)
                        {
                            fmpz_zero(padic_mat_entry(B, i, j));
                        }
                        else
                        {
                            fmpz_pow_ui(f, ctx->p, (v - w) - m);
                            fmpz_pow_ui(g, ctx->p, N - (v - w));
                            _padic_inv(t, t, ctx->p, N - (v - w));
                            fmpz_mul(padic_mat_entry(B, i, j), s, t);
                            fmpz_mod(padic_mat_entry(B, i, j), padic_mat_entry(B, i, j), g);
                            fmpz_mul(padic_mat_entry(B, i, j), padic_mat_entry(B, i, j), f);
                        }
                    }
            padic_mat_val(B) = m;
        }

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(s);
        fmpz_clear(t);
    }
}

