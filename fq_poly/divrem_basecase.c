/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

void _fq_poly_divrem_basecase(fq_struct *Q, fq_struct *R, 
    const fq_struct *A, long lenA, const fq_struct *B, long lenB, 
    const fq_t invB, const fq_ctx_t ctx)
{
    long iQ, iR;

    if (R != A)
        _fq_poly_set(R, A, lenA);

    for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
    {
        if (fq_is_zero(R + iR))
            fq_zero(Q + iQ);
        else
        {
            fq_mul(Q + iQ, R + iR, invB, ctx);

            _fq_poly_scalar_submul_fq(R + iQ, B, lenB, Q + iQ, ctx);
        }
    }
}

void fq_poly_divrem_basecase(fq_poly_t Q, fq_poly_t R, 
                             const fq_poly_t A, const fq_poly_t B, 
                             const fq_ctx_t ctx)
{
    const long lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    fq_struct *q, *r;
    fq_t invB;

    if (lenA < lenB)
    {
        fq_poly_set(R, A);
        fq_poly_zero(Q);
        return;
    }

    fq_init(invB);
    fq_inv(invB, fq_poly_lead(B), ctx);

    if (Q == A || Q == B)
    {
        q = _fq_poly_init(lenQ);
    }
    else
    {
        fq_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }
    if (R == B)
    {
        r = _fq_poly_init(lenA);
    }
    else
    {
        fq_poly_fit_length(R, lenA);
        r = R->coeffs;
    }

    _fq_poly_divrem_basecase(q, r, A->coeffs, lenA,
                                   B->coeffs, lenB, invB, ctx);

    if (Q == A || Q == B)
    {
        _fq_poly_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fq_poly_set_length(Q, lenQ);
    }
    if (R == B)
    {
        _fq_poly_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenA;
        R->length = lenA;
    }
    _fq_poly_set_length(R, lenB - 1);
    _fq_poly_normalise(R);

    fq_clear(invB);
}

