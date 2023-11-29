/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_div(mp_ptr Q, mp_srcptr A, slong lenA,  mp_srcptr B, slong lenB, nmod_t mod)
{
    if (lenA == lenB)
    {
        Q[0] = nmod_div(A[lenA - 1], B[lenB - 1], mod);
    }
    else if (lenB == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, nmod_inv(B[0], mod), mod);
    }
    else
    {
        gr_ctx_t ctx;
        _gr_ctx_init_nmod(ctx, &mod);

        if (lenB <= 15 || lenA - lenB <= 15)
            GR_MUST_SUCCEED(_gr_poly_div_basecase(Q, A, lenA, B, lenB, ctx));
        else
            GR_MUST_SUCCEED(_gr_poly_div_newton(Q, A, lenA, B, lenB, ctx));
    }
}

void
nmod_poly_div(nmod_poly_t Q,
                 const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tQ;
    mp_ptr q;
    slong A_len, B_len;

    B_len = B->length;

    if (B_len == 0)
    {
        if (nmod_poly_modulus(B) == 1)
        {
            nmod_poly_set(Q, A);
            return;
        } else
        {
            flint_throw(FLINT_ERROR, "Exception (nmod_poly_divrem). Division by zero.\n");
        }
    }

    A_len = A->length;

    if (A_len < B_len)
    {
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2(tQ, A->mod.n, A_len - B_len + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, A_len - B_len + 1);
        q = Q->coeffs;
    }

    _nmod_poly_div(q, A->coeffs, A_len,
                            B->coeffs, B_len, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(tQ, Q);
        nmod_poly_clear(tQ);
    }

    Q->length = A_len - B_len + 1;
}
