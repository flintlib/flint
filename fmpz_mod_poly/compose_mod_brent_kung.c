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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

void
_fmpz_mod_poly_compose_mod_brent_kung(fmpz * res, const fmpz * poly1, long len1,
                              const fmpz * poly2, const fmpz * poly3, long len3, const fmpz_t p)
{
    fmpz_mat_t A, B, C;
    fmpz * t, * h, * tmp;
    long i, j, n, m;

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
    fmpz_one(A->rows[0]);
    _fmpz_vec_set(A->rows[1], poly2, n);
    tmp = _fmpz_vec_init(2 * n - 1);
    for (i = 2; i < m; i++)
    {
        _fmpz_mod_poly_mulmod(tmp, A->rows[i - 1], n, poly2, n, poly3, len3, p);
        _fmpz_vec_set(A->rows[i], tmp, n);
    }
    _fmpz_vec_clear(tmp, 2 * n - 1);

    fmpz_mat_mul(C, B, A);
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_mod(C->rows[i] + j, C->rows[i] + j, p);

    /* Evaluate block composition using the Horner scheme */
    _fmpz_vec_set(res, C->rows[m - 1], n);
    _fmpz_mod_poly_mulmod(h, A->rows[m - 1], n, poly2, n, poly3, len3, p);

    for (i = m - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_mulmod(t, res, n, h, n, poly3, len3, p);
        _fmpz_mod_poly_add(res, t, n, C->rows[i], n, p);
    }

    _fmpz_vec_clear(h, 2 * n - 1);
    _fmpz_vec_clear(t, 2 * n - 1);

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
}

void
fmpz_mod_poly_compose_mod_brent_kung(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                             const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long len3 = poly3->length;
    long len = len3 - 1;
    long vec_len = FLINT_MAX(len3 - 1, len2);

    fmpz * ptr2;
    fmpz_t inv3;

    if (len3 == 0)
    {
        printf("Exception: division by zero in "
               "fmpz_mod_poly_compose_mod_brent_kung\n");
        abort();
    }

    if (len1 >= len3)
    {
        printf("Exception: fmpz_mod_poly_compose_brent_kung: the degree of the"
               " first polynomial must be smaller than that of the modulus\n");
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

    if (res == poly3 || res == poly1)
    {
        fmpz_mod_poly_t tmp;
        fmpz_mod_poly_init(tmp, &res->p);
        fmpz_mod_poly_compose_mod_brent_kung(tmp, poly1, poly2, poly3);
        fmpz_mod_poly_swap(tmp, res);
        fmpz_mod_poly_clear(tmp);
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
        fmpz_invmod(inv3, poly3->coeffs + len, &res->p);
        _fmpz_mod_poly_rem(ptr2, poly2->coeffs, len2,
                                 poly3->coeffs, len3, inv3, &res->p);
        fmpz_clear(inv3);
    }

    fmpz_mod_poly_fit_length(res, len);
    _fmpz_mod_poly_compose_mod_brent_kung(res->coeffs,
             poly1->coeffs, len1, ptr2, poly3->coeffs, len3, &res->p);
    _fmpz_mod_poly_set_length(res, len);
    _fmpz_mod_poly_normalise(res);

    _fmpz_vec_clear(ptr2, vec_len);
}
