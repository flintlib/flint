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
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

int
_TEMPLATE(T, poly_divides) (
    TEMPLATE(T, struct) * Q,
    const TEMPLATE(T, struct) * A, slong lenA,
    const TEMPLATE(T, struct) * B, slong lenB,
    const TEMPLATE(T, t) invB,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * R;
    slong lenR = 0;

    R = _TEMPLATE(T, vec_init) (lenA, ctx);

    _TEMPLATE(T, poly_divrem) (Q, R, A, lenA, B, lenB, invB, ctx);

    TEMPLATE(CAP_T, VEC_NORM) (R, lenR, ctx);
    _TEMPLATE(T, vec_clear) (R, lenA, ctx);

    return (lenR == 0);
}

int
TEMPLATE(T, poly_divides) (TEMPLATE(T, poly_t) Q,
                           const TEMPLATE(T, poly_t) A,
                           const TEMPLATE(T, poly_t) B,
                           const TEMPLATE(T, ctx_t) ctx)
{
    if (TEMPLATE(T, poly_is_zero) (B, ctx))
    {
        TEMPLATE_PRINTF("Exception (%s_poly_divides).  B is zero.\n", T);
        flint_abort();
    }

    if (TEMPLATE(T, poly_is_zero) (A, ctx))
    {
        TEMPLATE(T, poly_zero) (Q, ctx);
        return 1;
    }
    if (TEMPLATE(T, poly_length) (A, ctx) < TEMPLATE(T, poly_length) (B, ctx))
    {
        return 0;
    }

    {
        const slong lenQ =
            TEMPLATE(T, poly_length) (A, ctx) - TEMPLATE(T, poly_length) (B,
                                                                          ctx)
            + 1;
        int ans;
        TEMPLATE(T, t) invB;

        TEMPLATE(T, init) (invB, ctx);
        TEMPLATE(T, inv) (invB, TEMPLATE(T, poly_lead) (B, ctx), ctx);

        if (Q == A || Q == B)
        {
            TEMPLATE(T, poly_t) T;

            TEMPLATE(T, poly_init2) (T, lenQ, ctx);
            ans = _TEMPLATE(T, poly_divides) (T->coeffs, A->coeffs, A->length,
                                              B->coeffs, B->length, invB, ctx);
            _TEMPLATE(T, poly_set_length) (T, lenQ, ctx);
            _TEMPLATE(T, poly_normalise) (T, ctx);
            TEMPLATE(T, poly_swap) (Q, T, ctx);
            TEMPLATE(T, poly_clear) (T, ctx);
        }
        else
        {
            TEMPLATE(T, poly_fit_length) (Q, lenQ, ctx);
            ans = _TEMPLATE(T, poly_divides) (Q->coeffs, A->coeffs, A->length,
                                              B->coeffs, B->length, invB, ctx);
            _TEMPLATE(T, poly_set_length) (Q, lenQ, ctx);
            _TEMPLATE(T, poly_normalise) (Q, ctx);
        }
        TEMPLATE(T, clear) (invB, ctx);

        return ans;
    }
}


#endif
