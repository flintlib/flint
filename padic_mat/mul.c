/*
    Copyright (C) 2011, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_mat.h"

void padic_mat_mul(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx)
{
    if (padic_mat_is_empty(C))
    {
        return;
    }

    if (padic_mat_is_zero(A) || padic_mat_is_zero(B))
    {
        padic_mat_zero(C);
    }
    else
    {
        fmpz_mat_mul(padic_mat(C), padic_mat(A), padic_mat(B));

        padic_mat_val(C) = padic_mat_val(A) + padic_mat_val(B);

        padic_mat_reduce(C, ctx);
    }
}

