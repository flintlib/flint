/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

#include "mpn_extras.h"

void
_nmod_poly_divrem(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA, 
                                  mp_srcptr B, slong lenB, nmod_t mod)
{
    TMP_INIT;
    
    if (lenA == lenB)
        _nmod_poly_divrem_q0(Q, R, A, B, lenB, mod);
    else if (lenA == lenB + 1)
        _nmod_poly_divrem_q1(Q, R, A, lenA, B, lenB, mod);
    else if (lenB < 15)
    {
        mp_ptr W;
        
        TMP_START;
        W = TMP_ALLOC(NMOD_DIVREM_BC_ITCH(lenA, lenB, mod)*sizeof(mp_limb_t));
        _nmod_poly_divrem_basecase(Q, R, W, A, lenA, B, lenB, mod);
        TMP_END;
    }
    else if (lenB < 6000)
        _nmod_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, mod);
    else
        _nmod_poly_divrem_newton(Q, R, A, lenA, B, lenB, mod);
}

void nmod_polydr_divrem(nmod_polydr_t Q, nmod_polydr_t R,
            const nmod_polydr_t A, const nmod_polydr_t B, const nmod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length;
    nmod_polydr_t tQ, tR;
    mp_ptr q, r;
    
    if (lenB == 0)
    {
        flint_printf("Exception (nmod_poly_divrem). Division by zero.");
        flint_abort();
    }

    if (lenA < lenB)
    {
        nmod_polydr_set(R, A, ctx);
        nmod_polydr_zero(Q, ctx);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_polydr_init2(tQ, lenA - lenB + 1, ctx);
        q = tQ->coeffs;
    }
    else
    {
        nmod_polydr_fit_length(Q, lenA - lenB + 1, ctx);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_polydr_init2(tR, lenB - 1, ctx);
        r = tR->coeffs;
    }
    else
    {
        nmod_polydr_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    _nmod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, ctx->mod);

    if (Q == A || Q == B)
    {
        nmod_polydr_swap(Q, tQ, ctx);
        nmod_polydr_clear(tQ, ctx);
    }
    if (R == A || R == B)
    {
        nmod_polydr_swap(R, tR, ctx);
        nmod_polydr_clear(tR, ctx);
    }
        
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_polydr_normalise(R);
}

void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R,
                      const nmod_poly_t A, const nmod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    nmod_poly_t tQ, tR;
    mp_ptr q, r;
    
    if (lenB == 0)
    {
        flint_printf("Exception (nmod_poly_divrem). Division by zero.");
        flint_abort();
    }

    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2_preinv(tQ, A->mod.n, A->mod.ninv, lenA - lenB + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_poly_init2_preinv(tR, B->mod.n, B->mod.ninv, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, tQ);
        nmod_poly_clear(tQ);
    }
    if (R == A || R == B)
    {
        nmod_poly_swap(R, tR);
        nmod_poly_clear(tR);
    }
        
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}

