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

   Copyright (C) 2011 Fredrik Johansson
   Copyright (C) 2012 Lina Kulakova
   Copyright (C) 2013 Martin Lee

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

void
_fmpz_mod_poly_reduce_matrix_mod_poly (fmpz_mat_t A, const fmpz_mat_t B,
                                   const fmpz_mod_poly_t f)
{
    fmpz * tmp1, *tmp2;
    slong n = f->length - 1;
    slong i, m = n_sqrt(n) + 1;

    fmpz_t invf;
    fmpz_init(invf);
    fmpz_invmod(invf, f->coeffs + n, &f->p);

    fmpz_mat_init(A, m, n);

    fmpz_one(A->rows[0]);
    tmp1 = _fmpz_vec_init(2 * (B->c) - n);
    tmp2 = tmp1 + (B->c - n);
    for (i= 1; i < m; i++)
    {
        _fmpz_mod_poly_divrem(tmp1, tmp2, B->rows[i], B->c, f->coeffs,
                           f->length, invf, &f->p);
        _fmpz_vec_set(A->rows[i], tmp2, n);
    }

    _fmpz_vec_clear(tmp1, 2 * (B->c) - n);
    fmpz_clear(invf);
}

void
_fmpz_mod_poly_precompute_matrix (fmpz_mat_t A, const fmpz * poly1,
                          const fmpz * poly2, slong len2, const fmpz * poly2inv,
                          slong len2inv, const fmpz_t p)
{
    /* Set rows of A to powers of poly1 */
    slong i, n, m;

    n = len2 - 1;

    m = n_sqrt(n) + 1;

    fmpz_one(A->rows[0]);
    _fmpz_vec_set(A->rows[1], poly1, n);
    for (i = 2; i < m; i++)
        _fmpz_mod_poly_mulmod_preinv(A->rows[i], A->rows[i - 1], n, poly1, n,
                                      poly2, len2, poly2inv, len2inv, p);
}

void
fmpz_mod_poly_precompute_matrix(fmpz_mat_t A, const fmpz_mod_poly_t poly1,
                    const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly2inv)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len = len2 - 1;
    slong vec_len = FLINT_MAX(len2 - 1, len1);
    slong m= n_sqrt(len) + 1;

    fmpz* ptr;
    fmpz_t inv2;

    if (len2 == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_precompute_matrix)."
                     "Division by zero.\n");
        abort();
    }

    if (A->r != m || A->c != len)
    {
        flint_printf("Exception (fmpz_mod_poly_precompute_matrix)."
                     " Wrong dimensions.\n");
        abort();
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
        fmpz_invmod(inv2, poly2->coeffs + len, &poly1->p);
        _fmpz_mod_poly_rem(ptr, poly1->coeffs, len1,
                                 poly2->coeffs, len2, inv2, &poly1->p);
        fmpz_clear(inv2);
    }

    _fmpz_mod_poly_precompute_matrix (A, ptr, poly2->coeffs, len2,
                                 poly2inv->coeffs, poly2inv->length, &poly1->p);

    _fmpz_vec_clear(ptr, vec_len);
}

void
_fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz * res,
         const fmpz * poly1, slong len1, const fmpz_mat_t A, const fmpz * poly3,
         slong len3, const fmpz * poly3inv, slong len3inv, const fmpz_t p)
{
    fmpz_mat_t B, C;
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
        _fmpz_mod_poly_evaluate_fmpz(res, poly1, len1, A->rows[1], p);
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
        for (j = 0; j < n; j++)
            fmpz_mod(C->rows[i] + j, C->rows[i] + j, p);

    /* Evaluate block composition using the Horner scheme */
    _fmpz_vec_set(res, C->rows[m - 1], n);
    _fmpz_mod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly3,
                                 len3, poly3inv, len3inv, p);

    for (i = m - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                                     poly3inv, len3inv, p);
        _fmpz_mod_poly_add(res, t, n, C->rows[i], n, p);
    }

    _fmpz_vec_clear(h, n);
    _fmpz_vec_clear(t, n);

    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
}

void
fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz_mod_poly_t res,
                    const fmpz_mod_poly_t poly1, const fmpz_mat_t A,
                    const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv)
{
    slong len1 = poly1->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    if (len3 == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv)."
                     "Division by zero\n");
        abort();
    }

    if (len1 >= len3)
    {
        flint_printf("Exception (fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv)."
               "The degree of the first polynomial must be smaller than that of the "
               " modulus\n");
        abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        fmpz_mod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        fmpz_mod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1 || res == poly3inv)
    {
        fmpz_mod_poly_t tmp;
        fmpz_mod_poly_init(tmp, &res->p);
        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(tmp, poly1, A,
                                                    poly3, poly3inv);
        fmpz_mod_poly_swap(tmp, res);
        fmpz_mod_poly_clear(tmp);
        return;
    }

    fmpz_mod_poly_fit_length(res, len);
    _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(res->coeffs,
             poly1->coeffs, len1, A, poly3->coeffs, len3,
             poly3inv->coeffs, poly3inv->length, &res->p);
    _fmpz_mod_poly_set_length(res, len);
    _fmpz_mod_poly_normalise(res);

}
