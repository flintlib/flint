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

int
_fq_poly_divides(fq_struct * Q,
                 const fq_struct * A, long lenA,
                 const fq_struct * B, long lenB, const fq_t invB,
                 const fq_ctx_t ctx)
{
    fq_struct *R;
    long lenR;

    R = _fq_vec_init(lenA, ctx);

    _fq_poly_divrem(Q, R, A, lenA, B, lenB, invB, ctx);

    FQ_VEC_NORM(R, lenR, ctx);
    _fq_vec_clear(R, lenA, ctx);

    return (lenR == 0);
}

int
fq_poly_divides(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B,
                const fq_ctx_t ctx)
{
    if (fq_poly_is_zero(B, ctx))
    {
        flint_printf("Exception (fq_poly_divides).  B is zero.\n");
        abort();
    }

    if (fq_poly_is_zero(A, ctx))
    {
        fq_poly_zero(Q, ctx);
        return 1;
    }
    if (fq_poly_length(A, ctx) < fq_poly_length(B, ctx))
    {
        return 0;
    }

    {
        const long lenQ = fq_poly_length(A, ctx) - fq_poly_length(B, ctx) + 1;
        int ans;
        fq_t invB;

        fq_init(invB, ctx);
        fq_inv(invB, fq_poly_lead(B, ctx), ctx);

        if (Q == A || Q == B)
        {
            fq_poly_t T;

            fq_poly_init2(T, lenQ, ctx);
            ans = _fq_poly_divides(T->coeffs, A->coeffs, A->length,
                                   B->coeffs, B->length, invB, ctx);
            _fq_poly_set_length(T, lenQ, ctx);
            _fq_poly_normalise(T, ctx);
            fq_poly_swap(Q, T, ctx);
            fq_poly_clear(T, ctx);
        }
        else
        {
            fq_poly_fit_length(Q, lenQ, ctx);
            ans = _fq_poly_divides(Q->coeffs, A->coeffs, A->length,
                                   B->coeffs, B->length, invB, ctx);
            _fq_poly_set_length(Q, lenQ, ctx);
            _fq_poly_normalise(Q, ctx);
        }
        fq_clear(invB, ctx);

        return ans;
    }
}
