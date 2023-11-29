/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly1, slong len1,
    const TEMPLATE(T, struct) * poly2,
    const TEMPLATE(T, struct) * poly3, slong len3,
    const TEMPLATE(T, struct) * poly3inv, slong len3inv,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, mat_t) A, B, C;
    TEMPLATE(T, struct) * t, *h, *tmp;
    slong i, n, m;

    n = len3 - 1;

    if (len3 == 1)
        return;

    if (len1 == 1)
    {
        TEMPLATE(T, set) (res, poly1, ctx);
        return;
    }

    if (len3 == 2)
    {
        _TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (res, poly1, len1, poly2,
                                                  ctx);
        return;
    }

    m = n_sqrt(n) + 1;

    TEMPLATE(T, mat_init) (A, m, n, ctx);
    TEMPLATE(T, mat_init) (B, m, m, ctx);
    TEMPLATE(T, mat_init) (C, m, n, ctx);

    h = _TEMPLATE(T, vec_init) (2 * n - 1, ctx);
    t = _TEMPLATE(T, vec_init) (2 * n - 1, ctx);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _TEMPLATE(T, vec_set) (B->rows[i], poly1 + i * m, m, ctx);

    _TEMPLATE(T, vec_set) (B->rows[i], poly1 + i * m, len1 % m, ctx);

    /* Set rows of A to powers of poly2 */
    TEMPLATE(T, one) (A->rows[0], ctx);
    _TEMPLATE(T, vec_set) (A->rows[1], poly2, n, ctx);
    tmp = _TEMPLATE(T, vec_init) (2 * n - 1, ctx);
    for (i = 2; i < m; i++)
    {
        _TEMPLATE(T, poly_mulmod_preinv) (tmp, A->rows[i - 1], n, poly2, n,
                                          poly3, len3, poly3inv, len3inv, ctx);
        _TEMPLATE(T, vec_set) (A->rows[i], tmp, n, ctx);
    }
    _TEMPLATE(T, vec_clear) (tmp, 2 * n - 1, ctx);

    TEMPLATE(T, mat_mul) (C, B, A, ctx);

    /* Evaluate block composition using the Horner scheme */
    _TEMPLATE(T, vec_set) (res, C->rows[m - 1], n, ctx);
    _TEMPLATE(T, poly_mulmod_preinv) (h, A->rows[m - 1], n, poly2, n, poly3,
                                      len3, poly3inv, len3inv, ctx);

    for (i = m - 2; i >= 0; i--)
    {
        _TEMPLATE(T, poly_mulmod_preinv) (t, res, n, h, n, poly3, len3,
                                          poly3inv, len3inv, ctx);
        _TEMPLATE(T, poly_add) (res, t, n, C->rows[i], n, ctx);
    }

    _TEMPLATE(T, vec_clear) (h, 2 * n - 1, ctx);
    _TEMPLATE(T, vec_clear) (t, 2 * n - 1, ctx);

    TEMPLATE(T, mat_clear) (A, ctx);
    TEMPLATE(T, mat_clear) (B, ctx);
    TEMPLATE(T, mat_clear) (C, ctx);
}

void
TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly1,
    const TEMPLATE(T, poly_t) poly2,
    const TEMPLATE(T, poly_t) poly3,
    const TEMPLATE(T, poly_t) poly3inv,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len3inv = poly3inv->length;
    slong len = len3 - 1;
    slong vec_len = FLINT_MAX(len3 - 1, len2);

    TEMPLATE(T, struct) * ptr2;
    TEMPLATE(T, t) inv3;

    if (len3 == 0)
    {
        flint_throw(FLINT_ERROR, "(%s): Division by zero\n", __func__);
    }

    if (len1 >= len3)
    {
        flint_throw(FLINT_ERROR, "(%s): The degree of the first polynomial must be smaller than that of the modulus\n", __func__);
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
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (tmp, poly1, poly2,
                                                         poly3, poly3inv, ctx);
        TEMPLATE(T, poly_swap) (tmp, res, ctx);
        TEMPLATE(T, poly_clear) (tmp, ctx);
        return;
    }

    ptr2 = _TEMPLATE(T, vec_init) (vec_len, ctx);

    if (len2 <= len)
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
    _TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (res->coeffs,
                                                      poly1->coeffs, len1,
                                                      ptr2, poly3->coeffs,
                                                      len3, poly3inv->coeffs,
                                                      len3inv, ctx);
    _TEMPLATE(T, poly_set_length) (res, len, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);

    _TEMPLATE(T, vec_clear) (ptr2, vec_len, ctx);
}


#endif
