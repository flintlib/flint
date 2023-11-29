/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "flint.h"
#include "ulong_extras.h"
void
_TEMPLATE(T, poly_compose_mod_horner) (
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * g,
    const TEMPLATE(T, struct) * h, slong lenh,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, len;
    TEMPLATE(T, struct) * t;

    if (lenh == 1)
        return;

    if (lenf == 1)
    {
        TEMPLATE(T, set) (res, f, ctx);
        return;
    }

    if (lenh == 2)
    {
        _TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (res, f, lenf, g, ctx);
        return;
    }

    len = lenh - 1;
    i = lenf - 1;
    t = _TEMPLATE(T, vec_init) (2 * lenh - 3, ctx);

    _TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (res, g, len, f + i, ctx);
    i--;
    if (i >= 0)
    {
        TEMPLATE(T, add) (res, res, f + i, ctx);
    }

    while (i > 0)
    {
        i--;
        _TEMPLATE(T, poly_mulmod) (t, res, len, g, len, h, lenh, ctx);
        _TEMPLATE(T, poly_add) (res, t, len, f + i, 1, ctx);
    }

    _TEMPLATE(T, vec_clear) (t, 2 * lenh - 3, ctx);
}

void
TEMPLATE(T, poly_compose_mod_horner) (TEMPLATE(T, poly_t) res,
                                      const TEMPLATE(T, poly_t) poly1,
                                      const TEMPLATE(T, poly_t) poly2,
                                      const TEMPLATE(T, poly_t) poly3,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) inv3;
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;
    slong vec_len = FLINT_MAX(len3 - 1, len2);

    TEMPLATE(T, struct) * ptr2;

    if (len3 == 0)
    {
        flint_throw(FLINT_ERROR, "(%s): Division by zero\n", __func__);
    }

    if (len1 == 0 || len3 == 1)
    {
        TEMPLATE(T, poly_zero) (res, ctx);
        return;
    }

    if (len1 == 1)
    {
        TEMPLATE(T, poly_set) (res, poly1, ctx);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        TEMPLATE(T, poly_t) tmp;
        TEMPLATE(T, poly_init) (tmp, ctx);
        TEMPLATE(T, poly_compose_mod_horner) (tmp, poly1, poly2, poly3, ctx);
        TEMPLATE(T, poly_swap) (tmp, res, ctx);
        TEMPLATE(T, poly_clear) (tmp, ctx);
        return;
    }

    ptr2 = _TEMPLATE(T, vec_init) (vec_len, ctx);

    if (len2 <= len3 - 1)
    {
        _TEMPLATE(T, vec_set) (ptr2, poly2->coeffs, len2, ctx);
        _TEMPLATE(T, vec_zero) (ptr2 + len2, vec_len - len2, ctx);
    }
    else
    {
        TEMPLATE(T, init) (inv3, ctx);
        TEMPLATE(T, inv) (inv3, poly3->coeffs + len, ctx);
        _TEMPLATE(T, poly_rem) (ptr2, poly2->coeffs, len2,
                                poly3->coeffs, len3, inv3, ctx);
        TEMPLATE(T, clear) (inv3, ctx);
    }

    TEMPLATE(T, poly_fit_length) (res, len, ctx);
    _TEMPLATE(T, poly_compose_mod_horner) (res->coeffs,
                                           poly1->coeffs, len1,
                                           ptr2, poly3->coeffs, len3, ctx);
    _TEMPLATE(T, poly_set_length) (res, len, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);

    _TEMPLATE(T, vec_clear) (ptr2, vec_len, ctx);
}


#endif
