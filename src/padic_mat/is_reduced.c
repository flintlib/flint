/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic.h"
#include "padic_mat.h"

int
padic_mat_is_reduced(const padic_mat_t A, const padic_ctx_t ctx)
{
    if (padic_mat_is_empty(A))
    {
        return 1;
    }
    else if (fmpz_mat_is_zero(padic_mat(A)))
    {
        return (padic_mat_val(A) == 0);
    }
    else if (padic_mat_is_canonical(A, ctx))
    {
        const slong v = padic_mat_val(A);
        const slong N = padic_mat_prec(A);

        if (v >= N)
        {
            return 0;
        }
        else
        {
            slong i, j;
            fmpz_t pN;
            int reduced = 1;
            int alloc = _padic_ctx_pow_ui(pN, N - v, ctx);

            for (i = 0; (i < padic_mat_nrows(A)) && reduced; i++)
                for (j = 0; (j < padic_mat_ncols(A)) && reduced; j++)
                    reduced = (fmpz_cmp(padic_mat_entry(A, i, j), pN) < 0);

            if (alloc)
                fmpz_clear(pN);

            return reduced;
        }
    }
    else
    {
        return 0;
    }
}
