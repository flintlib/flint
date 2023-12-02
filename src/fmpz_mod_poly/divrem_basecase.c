/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_divrem_basecase(fmpz *Q, fmpz *R,
    const fmpz *A, slong lenA, const fmpz *B, slong lenB,
    const fmpz_t invB, const fmpz_mod_ctx_t ctx)
{
    slong iQ, iR;
    fmpz * W;
    TMP_INIT;

	TMP_START;

    if (R != A)
    {
        /* cannot use R as it might not have enough space */
        FMPZ_VEC_TMP_INIT(W, lenA);
        _fmpz_vec_set(W, A, lenA);
    }
    else
    {
       W = R;
    }

    for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
    {
        if (fmpz_is_zero(W + iR))
        {
            fmpz_zero(Q + iQ);
		}
        else
        {
            fmpz_mul(Q + iQ, W + iR, invB);
            fmpz_mod_set_fmpz(Q + iQ, Q + iQ, ctx);
            _fmpz_vec_scalar_submul_fmpz(W + iQ, B, lenB, Q + iQ);
        }

        if (iQ > 0)
            fmpz_mod_set_fmpz(W + iR - 1, W + iR - 1, ctx);
    }

	_fmpz_mod_vec_set_fmpz_vec(W, W, lenB - 1, ctx);

    if (R != A)
    {
        _fmpz_vec_swap(R, W, lenB - 1);
        FMPZ_VEC_TMP_CLEAR(W, lenA);
    }

	TMP_END;
}

void fmpz_mod_poly_divrem_basecase(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    fmpz *q, *r;
    fmpz_t invB;

    if (lenB == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A, ctx);
            fmpz_mod_poly_zero(R, ctx);
            return;
        }
        else
        {
            flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_divrem_basecase). Division by zero.\n");
        }
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A, ctx);
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    fmpz_init(invB);
    fmpz_invmod(invB, B->coeffs + (lenB - 1), fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }
    if (R == B)
    {
        r = _fmpz_vec_init(lenB - 1);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    _fmpz_mod_poly_divrem_basecase(q, r, A->coeffs, lenA, B->coeffs, lenB, invB, ctx);

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fmpz_mod_poly_set_length(Q, lenQ);
    }
    if (R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenB - 1;
        R->length = lenB - 1;
    }
    else
    {
      _fmpz_mod_poly_set_length(R, lenB - 1);
    }

    _fmpz_mod_poly_normalise(R);

    fmpz_clear(invB);
}

