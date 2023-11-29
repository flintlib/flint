/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013, 2014 Martin Lee
    Copyright (C) 2020 William Hart

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
_fmpz_mod_poly_reduce_matrix_mod_poly(fmpz_mat_t A, const fmpz_mat_t B,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    fmpz * tmp1, *tmp2;
    slong n = f->length - 1;
    slong i, m = n_sqrt(n) + 1;

    fmpz_t invf;
    fmpz_init(invf);
    fmpz_invmod(invf, f->coeffs + n, fmpz_mod_ctx_modulus(ctx));

    fmpz_mat_init(A, m, n);

    tmp1 = _fmpz_vec_init(2 * (B->c) - n);
    tmp2 = tmp1 + (B->c - n);

    fmpz_one(A->rows[0]);

    for (i= 1; i < m; i++)
    {
        _fmpz_mod_poly_divrem(tmp1, tmp2, B->rows[i], B->c, f->coeffs,
                                   f->length, invf, ctx);
        _fmpz_vec_set(A->rows[i], tmp2, n);
    }

    _fmpz_vec_clear(tmp1, 2*(B->c) - n);
    fmpz_clear(invf);
}

void
_fmpz_mod_poly_precompute_matrix(fmpz_mat_t A, const fmpz * poly1,
                          const fmpz * poly2, slong len2, const fmpz * poly2inv,
                          slong len2inv, const fmpz_mod_ctx_t ctx)
{
    /* Set rows of A to powers of poly1 */
    slong n, m;

    n = len2 - 1;

    m = n_sqrt(n) + 1;

    _fmpz_mod_poly_powers_mod_preinv_naive(A->rows, poly1, n, m,
                                            poly2, len2, poly2inv, len2inv, ctx);
}

void
fmpz_mod_poly_precompute_matrix(fmpz_mat_t A, const fmpz_mod_poly_t poly1,
                  const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly2inv,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len = len2 - 1;
    slong vec_len = FLINT_MAX(len2 - 1, len1);
    slong m = n_sqrt(len) + 1;

    fmpz * ptr;
    fmpz_t inv2;

    if (len2 == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mod_poly_precompute_matrix): Division by zero.\n");
    }

    if (A->r != m || A->c != len)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mod_poly_precompute_matrix): Wrong dimensions.\n");
    }

    if (len2 == 1)
    {
        fmpz_mat_zero(A);

        return;
    }

    ptr = _fmpz_vec_init(vec_len);

    if (len1 <= len)
    {
        _fmpz_vec_set(ptr, poly1->coeffs, len1);
        _fmpz_vec_zero(ptr + len1, vec_len - len1);
    }
    else
    {
        fmpz_init(inv2);
        fmpz_invmod(inv2, poly2->coeffs + len, fmpz_mod_ctx_modulus(ctx));
        _fmpz_mod_poly_rem(ptr, poly1->coeffs, len1,
                         poly2->coeffs, len2, inv2, ctx);
        fmpz_clear(inv2);
    }

    _fmpz_mod_poly_precompute_matrix (A, ptr, poly2->coeffs, len2,
                poly2inv->coeffs, poly2inv->length, ctx);

    _fmpz_vec_clear(ptr, vec_len);
}

void
_fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz * res,
         const fmpz * poly1, slong len1, const fmpz_mat_t A, const fmpz * poly3,
         slong len3, const fmpz * poly3inv, slong len3inv, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_t B, C;
    fmpz * t, * h;
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
        _fmpz_mod_poly_evaluate_fmpz(res, poly1, len1, A->rows[1], ctx);
        return;
    }

    m = n_sqrt(n) + 1;

    fmpz_mat_init(B, m, m);
    fmpz_mat_init(C, m, n);

    h = _fmpz_vec_init(n);
    t = _fmpz_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _fmpz_vec_set(B->rows[i], poly1 + i * m, m);

    _fmpz_vec_set(B->rows[i], poly1 + i * m, len1 % m);

    fmpz_mat_mul(C, B, A);
    for (i = 0; i < m; i++)
        _fmpz_mod_vec_set_fmpz_vec(C->rows[i], C->rows[i], n, ctx);

    /* Evaluate block composition using the Horner scheme */
    _fmpz_vec_set(res, C->rows[m - 1], n);
    _fmpz_mod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly3,
                                 len3, poly3inv, len3inv, ctx);

    for (i = m - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                                     poly3inv, len3inv, ctx);
        _fmpz_mod_poly_add(res, t, n, C->rows[i], n, ctx);
    }

    _fmpz_vec_clear(h, n);
    _fmpz_vec_clear(t, n);

    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
}

void fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz_mod_poly_t res,
                  const fmpz_mod_poly_t poly1, const fmpz_mat_t A,
                  const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    if (len3 == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv): "
                     "Division by zero\n");
    }

    if (len1 >= len3)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv): "
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

    if (res == poly3 || res == poly1 || res == poly3inv)
    {
        fmpz_mod_poly_t tmp;

        fmpz_mod_poly_init(tmp, ctx);

        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(tmp, poly1, A,
                                                         poly3, poly3inv, ctx);
        fmpz_mod_poly_swap(tmp, res, ctx);
        fmpz_mod_poly_clear(tmp, ctx);

        return;
    }

    fmpz_mod_poly_fit_length(res, len, ctx);

    _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(res->coeffs,
             poly1->coeffs, len1, A, poly3->coeffs, len3,
             poly3inv->coeffs, poly3inv->length, ctx);

    _fmpz_mod_poly_set_length(res, len);
    _fmpz_mod_poly_normalise(res);
}
