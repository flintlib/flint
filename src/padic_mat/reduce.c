/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_mat.h"

void _padic_mat_reduce(padic_mat_t mat, const padic_ctx_t ctx)
{
    if (!padic_mat_is_empty(mat) && !padic_mat_is_zero(mat))
    {
        if (mat->val >= padic_mat_prec(mat))
        {
            padic_mat_zero(mat);
        }
        else
        {
            slong i;
            fmpz_t pow;

            fmpz_init(pow);
            fmpz_pow_ui(pow, ctx->p, padic_mat_prec(mat) - mat->val);
            for (i = 0; i < padic_mat(mat)->r * padic_mat(mat)->c; i++)
            {
                fmpz_mod(padic_mat(mat)->entries + i, 
                         padic_mat(mat)->entries + i, pow);
            }
            fmpz_clear(pow);

            if (padic_mat_is_zero(mat))
            {
                mat->val = 0;
            }
        }
    }
}

void padic_mat_reduce(padic_mat_t mat, const padic_ctx_t ctx)
{
    _padic_mat_canonicalise(mat, ctx);
    _padic_mat_reduce(mat, ctx);
}

