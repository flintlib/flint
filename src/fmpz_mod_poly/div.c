/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart
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

void _fmpz_mod_poly_div(fmpz * Q,
    const fmpz *A, slong lenA, const fmpz *B, slong lenB,
    const fmpz_t invB, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);

    if (lenB <= 15 || lenA - lenB <= 15)
        GR_MUST_SUCCEED(_gr_poly_div_basecase_preinv1(Q, A, lenA, B, lenB, invB, gr_ctx));
    else
        GR_MUST_SUCCEED(_gr_poly_div_newton(Q, A, lenA, B, lenB, gr_ctx));  /* todo: pass inverse */
}

void fmpz_mod_poly_div(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    fmpz *q;
    fmpz_t invB;

    if (lenB == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A, ctx);
            return;
        }
        else
        {
            flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_div_basecase). Division by zero.\n");
        }
    }

    if (lenA < lenB)
    {
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

    _fmpz_mod_poly_div(q, A->coeffs, lenA, B->coeffs, lenB, invB, ctx);

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

    fmpz_clear(invB);
}

