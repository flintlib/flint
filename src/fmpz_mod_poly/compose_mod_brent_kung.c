/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_compose_mod_brent_kung(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, const fmpz * poly3, slong len3,
                                                                    const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_t A, B, C;
    fmpz * t, * h, * tmp;
    slong i, n, m;

    n = len3 - 1;

    if (len3 == 1)
        return;

    if (len1 == 1)
    {
        fmpz_set(res, poly1);
        return;
    }

    if (len3 == 2)
    {
        _fmpz_mod_poly_evaluate_fmpz(res, poly1, len1, poly2, ctx);
        return;
    }

    m = n_sqrt(n) + 1;

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, m, m);
    fmpz_mat_init(C, m, n);

    h = _fmpz_vec_init(2 * n - 1);
    t = _fmpz_vec_init(2 * n - 1);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _fmpz_vec_set(B->rows[i], poly1 + i * m, m);

    _fmpz_vec_set(B->rows[i], poly1 + i * m, len1 % m);

    /* Set rows of A to powers of poly2 */
    fmpz_one(A->rows[0]);
    _fmpz_vec_set(A->rows[1], poly2, n);
    tmp = _fmpz_vec_init(2 * n - 1);
    for (i = 2; i < m; i++)
    {
        _fmpz_mod_poly_mulmod(tmp, A->rows[i - 1], n, poly2, n, poly3, len3, ctx);
        _fmpz_vec_set(A->rows[i], tmp, n);
    }
    _fmpz_vec_clear(tmp, 2 * n - 1);

    fmpz_mat_mul(C, B, A);
    for (i = 0; i < m; i++)
        _fmpz_mod_vec_set_fmpz_vec(C->rows[i], C->rows[i], n, ctx);

    /* Evaluate block composition using the Horner scheme */
    _fmpz_vec_set(res, C->rows[m - 1], n);
    _fmpz_mod_poly_mulmod(h, A->rows[m - 1], n, poly2, n, poly3, len3, ctx);

    for (i = m - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_mulmod(t, res, n, h, n, poly3, len3, ctx);
        _fmpz_mod_poly_add(res, t, n, C->rows[i], n, ctx);
    }

    _fmpz_vec_clear(h, 2 * n - 1);
    _fmpz_vec_clear(t, 2 * n - 1);

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
}

void fmpz_mod_poly_compose_mod_brent_kung(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                         const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;
    slong vec_len = FLINT_MAX(len3 - 1, len2);

    fmpz * ptr2;
    fmpz_t inv3;

    if (len3 == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mod_poly_compose_mod_brent_kung): Division by zero\n");
    }

    if (len1 >= len3)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mod_poly_compose_mod_brent_kung): "
                "The degree of the first polynomial must be smaller than that of the modulus\n");
    }

    if (len1 == 0 || len3 == 1)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    if (len1 == 1)
    {
        fmpz_mod_poly_set(res, poly1, ctx);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        fmpz_mod_poly_t tmp;
        fmpz_mod_poly_init(tmp, ctx);
        fmpz_mod_poly_compose_mod_brent_kung(tmp, poly1, poly2, poly3, ctx);
        fmpz_mod_poly_swap(tmp, res, ctx);
        fmpz_mod_poly_clear(tmp, ctx);
        return;
    }

    ptr2 = _fmpz_vec_init(vec_len);

    if (len2 <= len)
    {
        _fmpz_vec_set(ptr2, poly2->coeffs, len2);
        _fmpz_vec_zero(ptr2 + len2, vec_len - len2);
    }
    else
    {
        fmpz_init(inv3);
        fmpz_invmod(inv3, poly3->coeffs + len, fmpz_mod_ctx_modulus(ctx));
        _fmpz_mod_poly_rem(ptr2, poly2->coeffs, len2,
                         poly3->coeffs, len3, inv3, ctx);
        fmpz_clear(inv3);
    }

    fmpz_mod_poly_fit_length(res, len, ctx);
    _fmpz_mod_poly_compose_mod_brent_kung(res->coeffs, poly1->coeffs, len1,
                         ptr2, poly3->coeffs, len3, ctx);
    _fmpz_mod_poly_set_length(res, len);
    _fmpz_mod_poly_normalise(res);

    _fmpz_vec_clear(ptr2, vec_len);
}
