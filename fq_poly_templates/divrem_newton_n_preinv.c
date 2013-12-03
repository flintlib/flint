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

    Copyright (C) 2011 William Hart
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_divrem_newton_n_preinv) (
    TEMPLATE(T, struct) * Q,
    TEMPLATE(T, struct) * R,
    const TEMPLATE(T, struct) * A, slong lenA,
    const TEMPLATE(T, struct) * B, slong lenB,
    const TEMPLATE(T, struct) * Binv, slong lenBinv,
    const TEMPLATE(T, ctx_t) ctx)
{
    const slong lenQ = lenA - lenB + 1;

    _TEMPLATE(T, poly_div_newton_n_preinv) (Q, A, lenA, B, lenB, Binv, lenBinv,
                                            ctx);

    if (lenB > 1)
    {
        if (lenQ >= lenB - 1)
            _TEMPLATE(T, poly_mullow) (R, Q, lenQ, B, lenB - 1, lenB - 1, ctx);
        else
            _TEMPLATE(T, poly_mullow) (R, B, lenB - 1, Q, lenQ, lenB - 1, ctx);

        _TEMPLATE(T, vec_sub) (R, A, R, lenB - 1, ctx);
    }
}

void
TEMPLATE(T, poly_divrem_newton_n_preinv) (TEMPLATE(T, poly_t) Q,
                                          TEMPLATE(T, poly_t) R,
                                          const TEMPLATE(T, poly_t) A,
                                          const TEMPLATE(T, poly_t) B,
                                          const TEMPLATE(T, poly_t) Binv,
                                          const TEMPLATE(T, ctx_t) ctx)
{
    const slong lenA = A->length, lenB = B->length, lenBinv = Binv->length;
    TEMPLATE(T, struct) * q, *r;

    if (lenB == 0)
    {
        TEMPLATE_PRINTF
            ("Exception (%s_poly_divrem_newton_n_preinv). Division by zero.\n",
             T);
        abort();
    }

    if (lenA < lenB)
    {
        TEMPLATE(T, poly_set) (R, A, ctx);
        TEMPLATE(T, poly_zero) (Q, ctx);
        return;
    }

    if (lenA > 2 * lenB - 2)
    {
        TEMPLATE_PRINTF("Exception (%s_poly_divrem_newton_n_preinv).\n", T);
    }

    if (Q == A || Q == B || Q == Binv)
    {
        q = _TEMPLATE(T, vec_init) (lenA - lenB + 1, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (Q, lenA - lenB + 1, ctx);
        q = Q->coeffs;
    }
    if (R == A || R == B || R == Binv)
    {
        r = _TEMPLATE(T, vec_init) (lenB - 1, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (R, lenB - 1, ctx);
        r = R->coeffs;
    }

    _TEMPLATE(T, poly_divrem_newton_n_preinv) (q, r, A->coeffs, lenA,
                                               B->coeffs, lenB, Binv->coeffs,
                                               lenBinv, ctx);

    if (Q == A || Q == B || Q == Binv)
    {
        _TEMPLATE(T, vec_clear) (Q->coeffs, lenA - lenB + 1, ctx);
        Q->coeffs = q;
        Q->alloc = lenA - lenB + 1;
    }
    if (R == A || R == B || R == Binv)
    {
        _TEMPLATE(T, vec_clear) (R->coeffs, lenB - 1, ctx);
        R->coeffs = r;
        R->alloc = lenB - 1;
    }
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _TEMPLATE(T, poly_normalise) (R, ctx);
}


#endif
