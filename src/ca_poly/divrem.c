/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_divrem(ca_ptr Q, ca_ptr R, ca_srcptr A, slong lenA, ca_srcptr B, slong lenB, const ca_t invB, ca_ctx_t ctx)
{
    _ca_poly_divrem_basecase(Q, R, A, lenA, B, lenB, invB, ctx);
}

int
ca_poly_divrem(ca_poly_t Q, ca_poly_t R,
    const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    ca_ptr q, r;
    ca_t invB;

    if (lenB == 0)
        return 0;

    if (lenA < lenB)
    {
        if (ca_check_is_zero(B->coeffs + lenB - 1, ctx) == T_FALSE)
        {
            ca_poly_set(R, A, ctx);
            ca_poly_zero(Q, ctx);
            return 1;
        }
        else
        {
            return 0;
        }
    }

    ca_init(invB, ctx);
    ca_inv(invB, B->coeffs + lenB - 1, ctx);

    if (CA_IS_SPECIAL(invB))
    {
        ca_clear(invB, ctx);
        return 0;
    }

    if (Q == A || Q == B)
    {
        q = _ca_vec_init(lenQ, ctx);
    }
    else
    {
        ca_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    if (R == B)
    {
        r = _ca_vec_init(lenA, ctx);
    }
    else
    {
        ca_poly_fit_length(R, lenA, ctx);
        r = R->coeffs;
    }

    _ca_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, invB, ctx);

    if (Q == A || Q == B)
    {
        _ca_vec_clear(Q->coeffs, Q->alloc, ctx);
        Q->coeffs = q;
        Q->alloc = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _ca_poly_set_length(Q, lenQ, ctx);
    }

    if (R == B)
    {
        _ca_vec_clear(R->coeffs, R->alloc, ctx);
        R->coeffs = r;
        R->alloc = lenA;
        R->length = lenA;
    }

    _ca_poly_set_length(R, lenB - 1, ctx);
    _ca_poly_normalise(R, ctx);

    ca_clear(invB, ctx);

    return 1;
}

int
ca_poly_div(ca_poly_t Q, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
{
    ca_poly_t R;
    int success;
    ca_poly_init(R, ctx);
    success = ca_poly_divrem(Q, R, A, B, ctx);
    ca_poly_clear(R, ctx);
    return success;
}

int
ca_poly_rem(ca_poly_t R, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
{
    ca_poly_t Q;
    int success;
    ca_poly_init(Q, ctx);
    success = ca_poly_divrem(Q, R, A, B, ctx);
    ca_poly_clear(Q, ctx);
    return success;
}
