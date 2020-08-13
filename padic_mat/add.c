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

/*
    Assumptions:

    o  That the matrix dimensions be compatible.
    o  That the matrices be non-empty.
    o  That ord_p(A) >= ord_p(B).
 */

void _padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                   const padic_ctx_t ctx)
{
    if (padic_mat_is_zero(A))
    {
        padic_mat_set(C, B, ctx);
        return;
    }
    if (padic_mat_is_zero(B))
    {
        padic_mat_set(C, A, ctx);
        return;
    }
    if (padic_mat_val(B) >= padic_mat_prec(C))
    {
        padic_mat_zero(C);
        return;
    }
    

    if (padic_mat_val(A) == padic_mat_val(B))
    {
        fmpz_mat_add(padic_mat(C), padic_mat(A), padic_mat(B));
        padic_mat_val(C) = padic_mat_val(B);
        _padic_mat_canonicalise(C, ctx);
    }
    else  /* padic_mat_val(A) > padic_mat_val(B) */
    {
        fmpz_t x;

        fmpz_init(x);
        fmpz_pow_ui(x, ctx->p, padic_mat_val(A) - padic_mat_val(B));

        if (C == B)
        {
            fmpz_mat_scalar_addmul_fmpz(padic_mat(C), padic_mat(A), x);
        }
        else if (C == A)
        {
            fmpz_mat_scalar_mul_fmpz(padic_mat(C), padic_mat(A), x);
            fmpz_mat_add(padic_mat(C), padic_mat(B), padic_mat(C));
            padic_mat_val(C) = padic_mat_val(B);
        }
        else
        {
            fmpz_mat_set(padic_mat(C), padic_mat(B));
            fmpz_mat_scalar_addmul_fmpz(padic_mat(C), padic_mat(A), x);
            padic_mat_val(C) = padic_mat_val(B);
        }

        fmpz_clear(x);
    }

    /* Reduction */
    {
        fmpz_t pow;
        int alloc = _padic_ctx_pow_ui(pow, padic_mat_prec(C)- padic_mat_val(C), ctx);

        /* TODO: Improve, use input precision */
        _fmpz_vec_scalar_mod_fmpz(padic_mat(C)->entries, 
            padic_mat(C)->entries, padic_mat_nrows(C)*padic_mat_ncols(C), pow);

        if (fmpz_mat_is_zero(padic_mat(C)))
        {
            padic_mat_val(C) = 0;
        }


        if (alloc)
            fmpz_clear(pow);
    }
}

void padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx)
{
    if (padic_mat_is_empty(C))
    {
        return;
    }

    if (padic_mat_val(A) >= padic_mat_val(B))
    {
        _padic_mat_add(C, A, B, ctx);
    }
    else
    {
        _padic_mat_add(C, B, A, ctx);
    }
}

