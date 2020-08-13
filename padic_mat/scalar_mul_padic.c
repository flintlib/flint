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

void _padic_mat_scalar_mul_padic(padic_mat_t B, 
                                 const padic_mat_t A, const padic_t c, 
                                 const padic_ctx_t ctx)
{
    if (padic_mat_is_empty(B))
    {
        return;
    }

    if (padic_is_zero(c) || padic_mat_is_zero(A))
    {
        padic_mat_zero(B);
        return;
    }

    fmpz_mat_scalar_mul_fmpz(padic_mat(B), padic_mat(A), padic_unit(c));
    padic_mat_val(B) = padic_mat_val(A) + padic_val(c);
}

void padic_mat_scalar_mul_padic(padic_mat_t B, 
                                const padic_mat_t A, const padic_t c, 
                                const padic_ctx_t ctx)
{
    _padic_mat_scalar_mul_padic(B, A, c, ctx);
    _padic_mat_reduce(B, ctx);
}

