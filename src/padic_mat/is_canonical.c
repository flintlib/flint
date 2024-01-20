/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "padic_mat.h"

int
padic_mat_is_canonical(const padic_mat_t A, const padic_ctx_t ctx)
{
    if (fmpz_mat_is_zero(padic_mat(A)))
    {
        return (padic_mat_val(A) == 0);
    }
    else
    {
        slong i, j;
        int canonical = 0;

        for (i = 0; i < padic_mat(A)->r; i++)
            for (j = 0; j < padic_mat(A)->c; j++)
                if (!fmpz_divisible(padic_mat_entry(A, i, j), ctx->p))
                    canonical = 1;
        return canonical;
    }
}
