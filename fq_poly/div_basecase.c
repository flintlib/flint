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

    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

void _fq_poly_div_basecase(fq_struct *Q, fq_struct *R, 
                           const fq_struct *A, slong lenA,
                           const fq_struct *B, slong lenB, 
                           const fq_t invB, const fq_ctx_t ctx)
{
    const slong alloc = (R == NULL) ? lenA : 0;
    slong lenR = lenB - 1, iQ;

    if (alloc)
        R = _fq_vec_init(alloc, ctx);
    if (R != A)
        _fq_vec_set(R + lenR, A + lenR, lenA - lenR, ctx);

    for (iQ = lenA - lenB; iQ >= 0; iQ--)
    {
        if (fq_is_zero(R + lenA - 1, ctx))
        {
            fq_zero(Q + iQ, ctx);
        }
        else
        {
            fq_mul(Q + iQ, R + lenA - 1, invB, ctx);

            _fq_vec_scalar_submul_fq(R + lenA - lenR - 1, B, lenR, Q + iQ, ctx);
        }

        if (lenR - 1 >= iQ)
        {
            B++;
            lenR--;
        }

        lenA--;
    }

    if (alloc)
        _fq_vec_clear(R, alloc, ctx);
}

void fq_poly_div_basecase(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B,
                          const fq_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    fq_struct *q;
    fq_t invB;

    if (lenA < lenB)
    {
        fq_poly_zero(Q, ctx);
        return;
    }

    fq_init(invB, ctx);
    fq_inv(invB, B->coeffs + (lenB - 1), ctx);

    if (Q == A || Q == B)
    {
        q = _fq_vec_init(lenQ, ctx);
    }
    else
    {
        fq_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    _fq_poly_div_basecase(q, NULL, A->coeffs, lenA,
                          B->coeffs, lenB, invB, ctx);

    if (Q == A || Q == B)
    {
        _fq_vec_clear(Q->coeffs, Q->alloc, ctx);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fq_poly_set_length(Q, lenQ, ctx);
    }

    fq_clear(invB, ctx);
}

