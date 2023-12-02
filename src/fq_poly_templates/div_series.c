/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_div_series) (TEMPLATE(T, struct) * Q, const TEMPLATE(T, struct) * A, slong Alen,
    const TEMPLATE(T, struct) * B, slong Blen, slong n, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) u, d;

    TEMPLATE(T, init)(d, ctx);
    TEMPLATE(T, init)(u, ctx);

    if (!TEMPLATE(T, is_one)(B + 0, ctx))
       TEMPLATE(T, inv)(u, B + 0, ctx);
    else
       TEMPLATE(T, set_si)(u, 1, ctx);

    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (Blen == 1)
    {
        if (TEMPLATE(T, is_one)(B + 0, ctx))
            _TEMPLATE(T, vec_set)(Q, A, Alen, ctx);
        else
           _TEMPLATE3(T, poly_scalar_mul, T)(Q, A, Alen, u, ctx);

       _TEMPLATE(T, vec_zero)(Q + Alen, n - Alen, ctx);
    }
    else if (n < 16 || Blen < 10)
    {
        slong i, j;
        TEMPLATE(T, t) temp;

        TEMPLATE(T, init)(temp, ctx);

        if (TEMPLATE(T, is_one)(B + 0, ctx))
            TEMPLATE(T, set)(Q + 0, A + 0, ctx);
        else
           TEMPLATE(T, mul)(Q + 0, u, A + 0, ctx);

        for (i = 1; i < n; i++)
        {
            TEMPLATE(T, mul)(Q + i, B + 1, Q + i - 1, ctx);

            for (j = 2; j < FLINT_MIN(i + 1, Blen); j++)
            {
               TEMPLATE(T, mul)(temp, B + j, Q + i - j, ctx);
               TEMPLATE(T, add)(Q + i, Q + i, temp, ctx);
            }

            if (i < Alen)
               TEMPLATE(T, sub)(Q + i, A + i, Q + i, ctx);
            else
               TEMPLATE(T, neg)(Q + i, Q + i, ctx);

            if (!TEMPLATE(T, is_one)(B + 0, ctx))
               TEMPLATE(T, mul)(Q + i, Q + i, u, ctx);
        }

        TEMPLATE(T, clear)(temp, ctx);
    }
    else
    {
        TEMPLATE(T, struct) * B2, * Binv = _TEMPLATE(T, vec_init)(n, ctx);

        if (n > Blen)
        {
           B2 = _TEMPLATE(T, vec_init)(n, ctx);
           _TEMPLATE(T, vec_set)(B2, B, Blen, ctx);
        } else
           B2 = (TEMPLATE(T, struct) *) B;

        _TEMPLATE(T, poly_inv_series)(Binv, B2, n, u, ctx);
        _TEMPLATE(T, poly_mullow)(Q, Binv, n, A, Alen, n, ctx);

        _TEMPLATE(T, vec_clear)(Binv, n, ctx);
        if (n > Blen)
           _TEMPLATE(T, vec_clear)(B2, n, ctx);
    }

    TEMPLATE(T, clear)(d, ctx);
    TEMPLATE(T, clear)(u, ctx);
}

void TEMPLATE(T, poly_div_series)(TEMPLATE(T, poly_t) Q, const TEMPLATE(T, poly_t) A,
                                         const TEMPLATE(T, poly_t) B, slong n, const TEMPLATE(T, ctx_t) ctx)
{
    slong Alen = FLINT_MIN(A->length, n);
    slong Blen = FLINT_MIN(B->length, n);

    if (Blen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fq_poly_div_series). Division by zero.\n");
    }

    if (Alen == 0)
    {
        TEMPLATE(T, poly_zero)(Q, ctx);
        return;
    }

    if (Q == A || Q == B)
    {
        TEMPLATE(T, poly_t) t;
        TEMPLATE(T, poly_init2)(t, n, ctx);
        _TEMPLATE(T, poly_div_series)(t->coeffs, A->coeffs, Alen, B->coeffs, Blen, n, ctx);
        TEMPLATE(T, poly_swap)(Q, t, ctx);
        TEMPLATE(T, poly_clear)(t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length)(Q, n, ctx);
        _TEMPLATE(T, poly_div_series)(Q->coeffs, A->coeffs, Alen, B->coeffs, Blen, n, ctx);
    }

    _TEMPLATE(T, poly_set_length)(Q, n, ctx);
    _TEMPLATE(T, poly_normalise)(Q, ctx);
}

#endif
