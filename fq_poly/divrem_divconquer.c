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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

static void
__fq_poly_divrem_divconquer(fq_struct * Q, fq_struct * R,
                            const fq_struct * A, slong lenA,
                            const fq_struct * B, slong lenB,
                            const fq_t invB, const fq_ctx_t ctx)
{
    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 n1 - 1 by n1 division
         */

        const slong n1 = lenA - lenB + 1;
        const slong n2 = lenB - n1;

        const fq_struct *p1 = A + n2;
        const fq_struct *d1 = B + n2;
        const fq_struct *d2 = B;

        fq_struct *W = _fq_vec_init((2 * n1 - 1) + lenB - 1, ctx);

        fq_struct *d1q1 = R + n2;
        fq_struct *d2q1 = W + (2 * n1 - 1);

        _fq_poly_divrem_divconquer_recursive(Q, d1q1, W, p1, d1, n1,
                                             invB, ctx);

        /*
           Compute d2q1 = Q d2, of length lenB - 1
         */

        if (n1 >= n2)
            _fq_poly_mul(d2q1, Q, n1, d2, n2, ctx);
        else
            _fq_poly_mul(d2q1, d2, n2, Q, n1, ctx);

        /*
           Compute BQ = d1q1 * x^n1 + d2q1, of length lenB - 1; 
           then compute R = A - BQ
         */

        _fq_vec_swap(R, d2q1, n2, ctx);
        _fq_poly_add(R + n2, R + n2, n1 - 1, d2q1 + n2, n1 - 1, ctx);
        _fq_poly_sub(R, A, lenA, R, lenA, ctx);

        _fq_vec_clear(W, (2 * n1 - 1) + lenB - 1, ctx);
    }
    else                        /* lenA = 2 * lenB - 1 */
    {
        fq_struct *W = _fq_vec_init(lenA, ctx);

        _fq_poly_divrem_divconquer_recursive(Q, R, W, A, B, lenB, invB, ctx);

        _fq_poly_sub(R, A, lenB - 1, R, lenB - 1, ctx);

        _fq_vec_clear(W, lenA, ctx);
    }
}

void
_fq_poly_divrem_divconquer(fq_struct * Q, fq_struct * R,
                           const fq_struct * A, slong lenA,
                           const fq_struct * B, slong lenB,
                           const fq_t invB, const fq_ctx_t ctx)
{
    if (lenA <= 2 * lenB - 1)
    {
        __fq_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, invB, ctx);
    }
    else                        /* lenA > 2 * lenB - 1 */
    {
        slong shift, n = 2 * lenB - 1;
        fq_struct *QB, *W;

        _fq_vec_set(R, A, lenA, ctx);
        W = _fq_vec_init(2 * n, ctx);
        QB = W + n;

        while (lenA >= n)
        {
            shift = lenA - n;
            _fq_poly_divrem_divconquer_recursive(Q + shift, QB,
                                                 W, R + shift, B, lenB,
                                                 invB, ctx);
            _fq_poly_sub(R + shift, R + shift, n, QB, n, ctx);
            lenA -= lenB;
        }

        if (lenA >= lenB)
        {
            __fq_poly_divrem_divconquer(Q, W, R, lenA, B, lenB, invB, ctx);
            _fq_vec_swap(W, R, lenA, ctx);
        }

        _fq_vec_clear(W, 2 * n, ctx);
    }
}

void
fq_poly_divrem_divconquer(fq_poly_t Q, fq_poly_t R,
                          const fq_poly_t A, const fq_poly_t B,
                          const fq_ctx_t ctx)
{
    const slong lenA = A->length;
    const slong lenB = B->length;
    const slong lenQ = lenA - lenB + 1;

    fq_struct *q, *r;
    fq_t invB;

    if (lenA < lenB)
    {
        fq_poly_set(R, A, ctx);
        fq_poly_zero(Q, ctx);
        return;
    }

    fq_init(invB, ctx);
    fq_inv(invB, fq_poly_lead(B, ctx), ctx);

    if (Q == A || Q == B)
    {
        q = _fq_vec_init(lenQ, ctx);
    }
    else
    {
        fq_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = _fq_vec_init(lenA, ctx);
    }
    else
    {
        fq_poly_fit_length(R, lenA, ctx);
        r = R->coeffs;
    }

    _fq_poly_divrem_divconquer(q, r, A->coeffs, lenA,
                               B->coeffs, lenB, invB, ctx);

    if (Q == A || Q == B)
    {
        _fq_vec_clear(Q->coeffs, Q->alloc, ctx);
        Q->coeffs = q;
        Q->alloc = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fq_poly_set_length(Q, lenQ, ctx);
    }

    if (R == A || R == B)
    {
        _fq_vec_clear(R->coeffs, R->alloc, ctx);
        R->coeffs = r;
        R->alloc = lenA;
        R->length = lenA;
    }
    _fq_poly_set_length(R, lenB - 1, ctx);
    _fq_poly_normalise(R, ctx);

    fq_clear(invB, ctx);
}
