/*
    Copyright (C) 2011 Fredrik Johansson
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

#include "flint.h"
#include "ulong_extras.h"

void
_TEMPLATE(T, poly_reduce_matrix_mod_poly) (TEMPLATE(T, mat_t) A,
                                           const TEMPLATE(T, mat_t) B,
                                           const TEMPLATE(T, poly_t) f,
                                           const TEMPLATE(T, ctx_t) ctx)
{
    slong n = f->length - 1;
    slong i, m = n_sqrt(n) + 1;
    TEMPLATE(T, t) invf;

    TEMPLATE(T, mat_init) (A, m, n, ctx);

    TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (A, 0, 0), ctx);
    TEMPLATE(T, init) (invf, ctx);
    TEMPLATE(T, inv) (invf, f->coeffs + (f->length - 1), ctx);
    for (i = 1; i < m; i++)
        _TEMPLATE(T, poly_rem) (A->rows[i], B->rows[i], B->c, f->coeffs,
                                f->length, invf, ctx);
    TEMPLATE(T, clear) (invf, ctx);
}

void
_TEMPLATE(T, poly_precompute_matrix) (
    TEMPLATE(T, mat_t) A,
    const TEMPLATE(T, struct) * poly1,
    const TEMPLATE(T, struct) * poly2, slong len2,
    const TEMPLATE(T, struct) * poly2inv, slong len2inv,
    const TEMPLATE(T, ctx_t) ctx)
{
    /* Set rows of A to powers of poly1 */
    slong i, n, m;

    n = len2 - 1;

    m = n_sqrt(n) + 1;

    TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (A, 0, 0), ctx);
    _TEMPLATE(T, vec_set) (A->rows[1], poly1, n, ctx);
    for (i = 2; i < m; i++)
        _TEMPLATE(T, poly_mulmod_preinv) (A->rows[i], A->rows[i - 1],
                                          n, poly1, n, poly2, len2,
                                          poly2inv, len2inv, ctx);
}

void
TEMPLATE(T, poly_precompute_matrix) (TEMPLATE(T, mat_t) A,
                                     const TEMPLATE(T, poly_t) poly1,
                                     const TEMPLATE(T, poly_t) poly2,
                                     const TEMPLATE(T, poly_t) poly2inv,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len = len2 - 1;
    slong m = n_sqrt(len) + 1;

    TEMPLATE(T, struct) * ptr1;

    if (len2 == 0)
    {
        flint_throw(FLINT_ERROR, "(%s): Division by zero.\n", __func__);
    }

    if (A->r != m || A->c != len)
    {
        flint_throw(FLINT_ERROR, "(%s): Wrong dimensions.\n", __func__);
    }

    if (len2 == 1)
    {
        TEMPLATE(T, mat_zero) (A, ctx);
        return;
    }

    ptr1 = _TEMPLATE(T, vec_init) (len, ctx);

    if (len1 <= len)
    {
        _TEMPLATE(T, vec_set) (ptr1, poly1->coeffs, len1, ctx);
        _TEMPLATE(T, vec_zero) (ptr1 + len1, len - len1, ctx);
    }
    else
    {
        TEMPLATE(T, t) inv2;
        TEMPLATE(T, init) (inv2, ctx);
        TEMPLATE(T, inv) (inv2, poly2->coeffs + len2 - 1, ctx);
        _TEMPLATE(T, poly_rem) (ptr1, poly1->coeffs, len1,
                                poly2->coeffs, len2, inv2, ctx);
        TEMPLATE(T, clear) (inv2, ctx);
    }

    _TEMPLATE(T, poly_precompute_matrix) (A, ptr1, poly2->coeffs, len2,
                                          poly2inv->coeffs, poly2inv->length,
                                          ctx);

    _TEMPLATE(T, vec_clear) (ptr1, len, ctx);
}

void
_TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv) (
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly1, slong len1,
    const TEMPLATE(T, mat_t) A,
    const TEMPLATE(T, struct) * poly3, slong len3,
    const TEMPLATE(T, struct) * poly3inv, slong len3inv,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, mat_t) B, C;
    TEMPLATE(T, struct) * t, *h;
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
        _TEMPLATE3(T, poly_evaluate, T) (res, poly1, len1,
                                         TEMPLATE(T, mat_entry) (A, 1, 0),
                                         ctx);
        return;
    }

    m = n_sqrt(n) + 1;

    /* TODO check A */

    TEMPLATE(T, mat_init) (B, m, m, ctx);
    TEMPLATE(T, mat_init) (C, m, n, ctx);

    h = _TEMPLATE(T, vec_init) (n, ctx);
    t = _TEMPLATE(T, vec_init) (n, ctx);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _TEMPLATE(T, vec_set) (B->rows[i], poly1 + i * m, m, ctx);

    _TEMPLATE(T, vec_set) (B->rows[i], poly1 + i * m, len1 % m, ctx);

    TEMPLATE(T, mat_mul) (C, B, A, ctx);

    /* Evaluate block composition using the Horner scheme */
    _TEMPLATE(T, vec_set) (res, C->rows[m - 1], n, ctx);
    _TEMPLATE(T, poly_mulmod_preinv) (h, A->rows[m - 1], n, A->rows[1], n,
                                      poly3, len3, poly3inv, len3inv, ctx);

    for (i = m - 2; i >= 0; i--)
    {
        _TEMPLATE(T, poly_mulmod_preinv) (t, res, n, h, n, poly3, len3,
                                          poly3inv, len3inv, ctx);
        _TEMPLATE(T, poly_add) (res, t, n, C->rows[i], n, ctx);
    }

    _TEMPLATE(T, vec_clear) (h, n, ctx);
    _TEMPLATE(T, vec_clear) (t, n, ctx);

    TEMPLATE(T, mat_clear) (B, ctx);
    TEMPLATE(T, mat_clear) (C, ctx);
}

void
TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv) (
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly1, const TEMPLATE(T, mat_t) A,
    const TEMPLATE(T, poly_t) poly3, const TEMPLATE(T, poly_t) poly3inv,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong len1 = poly1->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    if (len3 == 0)
    {
        flint_throw(FLINT_ERROR, "(%s): Division by zero.\n", __func__);
    }

    if (len1 >= len3)
    {
        flint_throw(FLINT_ERROR, "(%s): The degree of the first polynomial must be smaller than that of the modulus.\n", __func__);
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

    if (res == poly3 || res == poly1 || res == poly3inv)
    {
        TEMPLATE(T, poly_t) tmp;
        TEMPLATE(T, poly_init) (tmp, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv) (tmp, poly1, A,
                                                                 poly3,
                                                                 poly3inv,
                                                                 ctx);
        TEMPLATE(T, poly_swap) (tmp, res, ctx);
        TEMPLATE(T, poly_clear) (tmp, ctx);
        return;
    }

    TEMPLATE(T, poly_fit_length) (res, len, ctx);
    _TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv) (res->coeffs,
                                                              poly1->coeffs,
                                                              len1, A,
                                                              poly3->coeffs,
                                                              len3,
                                                              poly3inv->coeffs,
                                                              poly3inv->length,
                                                              ctx);
    res->length = len;
    _TEMPLATE(T, poly_normalise) (res, ctx);

}
#endif
