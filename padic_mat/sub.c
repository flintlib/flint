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

void _padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                   const padic_ctx_t ctx)
{
    if (padic_mat_is_zero(A))
    {
        padic_mat_neg(C, B, ctx);
        return;
    }
    if (padic_mat_is_zero(B))
    {
        padic_mat_set(C, A, ctx);
        return;
    }
    if (FLINT_MIN(padic_mat_val(A), padic_mat_val(B)) >= padic_mat_prec(C))
    {
        padic_mat_zero(C);
        return;
    }

    if (padic_mat_val(A) == padic_mat_val(B))
    {
        fmpz_mat_sub(padic_mat(C), padic_mat(A), padic_mat(B));
        padic_mat_val(C) = padic_mat_val(A);
        _padic_mat_canonicalise(C, ctx);
    }
    else 
    {
        fmpz_t x;
        fmpz_init(x);

        if (padic_mat_val(A) < padic_mat_val(B))
        {
            fmpz_pow_ui(x, ctx->p, padic_mat_val(B) - padic_mat_val(A));

            if (C == A)
            {
                fmpz_mat_scalar_submul_fmpz(padic_mat(C), padic_mat(B), x);
            }
            else if (C == B)
            {
                fmpz_neg(x, x);
                fmpz_mat_scalar_mul_fmpz(padic_mat(C), padic_mat(B), x);
                fmpz_mat_add(padic_mat(C), padic_mat(A), padic_mat(C));
                padic_mat_val(C) = padic_mat_val(A);
            }
            else
            {
                fmpz_mat_set(padic_mat(C), padic_mat(A));
                fmpz_mat_scalar_submul_fmpz(padic_mat(C), padic_mat(B), x);
                padic_mat_val(C) = padic_mat_val(A);
            }
        }
        else  /* A->val > B->val */
        {
            fmpz_pow_ui(x, ctx->p, padic_mat_val(A) - padic_mat_val(B));

            if (C == B)
            {
                fmpz_mat_scalar_submul_fmpz(padic_mat(C), padic_mat(A), x);
                fmpz_mat_neg(padic_mat(C), padic_mat(C));
            }
            else 
            {
                fmpz_mat_scalar_mul_fmpz(padic_mat(C), padic_mat(A), x);
                fmpz_mat_sub(padic_mat(C), padic_mat(C), padic_mat(B));
                padic_mat_val(C) = padic_mat_val(B);
            }
        }
        fmpz_clear(x);
    }

}

void padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx)
{
    if (padic_mat_is_empty(C))
    {
        return;
    }

    _padic_mat_sub(C, A, B, ctx);
    _padic_mat_reduce(C, ctx);
}

