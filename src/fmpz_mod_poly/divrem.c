/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr_poly.h"

void _fmpz_mod_poly_divrem(fmpz *Q, fmpz *R,
    const fmpz *A, slong lenA, const fmpz *B, slong lenB,
    const fmpz_t invB, const fmpz_mod_ctx_t ctx)
{
    if (lenB <= 30 || lenA - lenB <= 10)
    {
        _fmpz_mod_poly_divrem_basecase(Q, R, A, lenA, B, lenB, invB, ctx);
    }
    else
    {
        gr_ctx_t gr_ctx;
        _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
        GR_MUST_SUCCEED(_gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, gr_ctx));
    }
}

void
fmpz_mod_poly_divrem(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                            const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length;
    const slong lenB = B->length;
    const slong lenQ = lenA - lenB + 1;

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
            flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_divrem). Division by zero.\n");
        }
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A, ctx);
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

	if (B->length < 8)
	{
        fmpz_mod_poly_divrem_basecase(Q, R, A, B, ctx);
        return;
    }

    fmpz_init(invB);
    fmpz_invmod(invB, fmpz_mod_poly_lead(B, ctx), fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = _fmpz_vec_init(lenB - 1);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    _fmpz_mod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, invB, ctx);

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

    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenB - 1;
        R->length = lenB - 1;
    }

    _fmpz_mod_poly_set_length(R, lenB - 1);
    _fmpz_mod_poly_normalise(R);

    fmpz_clear(invB);
}
