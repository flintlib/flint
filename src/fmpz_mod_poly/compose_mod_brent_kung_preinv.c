/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

void
_fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz * res, const fmpz * poly1,
                 slong len1, const fmpz * poly2, const fmpz * poly3, slong len3,
                 const fmpz * poly3inv, slong len3inv, const fmpz_t p)
{
    fmpz_mat_t A, B, C;
    fmpz * t, * h;
    slong i, j, n, m;

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
        _fmpz_mod_poly_evaluate_fmpz(res, poly1, len1, poly2, p);
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
    _fmpz_mod_poly_powers_mod_preinv_naive(A->rows, poly2, n,
                                         m, poly3, len3, poly3inv, len3inv, p);

    fmpz_mat_mul(C, B, A);
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_mod(C->rows[i] + j, C->rows[i] + j, p);

    /* Evaluate block composition using the Horner scheme */
    _fmpz_vec_set(res, C->rows[m - 1], n);
    _fmpz_mod_poly_mulmod_preinv(h, A->rows[m - 1], n, poly2, n, poly3, len3,
                                 poly3inv, len3inv, p);

    for (i = m - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                                     poly3inv, len3inv, p);
        _fmpz_mod_poly_add(res, t, n, C->rows[i], n, p);
    }

    _fmpz_vec_clear(h, 2 * n - 1);
    _fmpz_vec_clear(t, 2 * n - 1);

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
}

void
fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz_mod_poly_t res,
                    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                    const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    fmpz * ptr2;
    fmpz_t inv3;

    if (len3 == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_compose_mod_brent_kung preinv)."
                     "Division by zero\n");
        flint_abort();
    }

    if (len1 >= len3)
    {
        flint_printf("Exception (fmpz_mod_poly_compose_mod_brent_kung_preinv)."
               "The degree of the first polynomial must be smaller than that of the "
               " modulus\n");
        flint_abort();
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
        fmpz_mod_poly_compose_mod_brent_kung_preinv(tmp, poly1, poly2,
                                                    poly3, poly3inv, ctx);
        fmpz_mod_poly_swap(tmp, res, ctx);
        fmpz_mod_poly_clear(tmp, ctx);
        return;
    }

    ptr2 = _fmpz_vec_init(len);

    if (len2 <= len)
    {
        _fmpz_vec_set(ptr2, poly2->coeffs, len2);
        _fmpz_vec_zero(ptr2 + len2, len - len2);
    }
    else
    {
        fmpz_init(inv3);
        fmpz_invmod(inv3, poly3->coeffs + len, fmpz_mod_ctx_modulus(ctx));
        _fmpz_mod_poly_rem(ptr2, poly2->coeffs, len2,
                         poly3->coeffs, len3, inv3, fmpz_mod_ctx_modulus(ctx));
        fmpz_clear(inv3);
    }

    fmpz_mod_poly_fit_length(res, len, ctx);
    _fmpz_mod_poly_compose_mod_brent_kung_preinv(res->coeffs,
             poly1->coeffs, len1, ptr2, poly3->coeffs, len3,
             poly3inv->coeffs, poly3inv->length, fmpz_mod_ctx_modulus(ctx));
    _fmpz_mod_poly_set_length(res, len);
    _fmpz_mod_poly_normalise(res);

    _fmpz_vec_clear(ptr2, len);
}
