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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"
#include "fq_mat.h"

void
_fq_poly_compose_mod_brent_kung_preinv(fq_struct * res,
                                       const fq_struct * poly1, slong len1,
                                       const fq_struct * poly2,
                                       const fq_struct * poly3, slong len3,
                                       const fq_struct * poly3inv,
                                       slong len3inv, const fq_ctx_t ctx)
{
    fq_mat_t A, B, C;
    fq_struct *t, *h, *tmp;
    slong i, n, m;

    n = len3 - 1;

    if (len3 == 1)
        return;

    if (len1 == 1)
    {
        fq_set(res, poly1, ctx);
        return;
    }

    if (len3 == 2)
    {
        _fq_poly_evaluate_fq(res, poly1, len1, poly2, ctx);
        return;
    }

    m = n_sqrt(n) + 1;

    fq_mat_init(A, m, n, ctx);
    fq_mat_init(B, m, m, ctx);
    fq_mat_init(C, m, n, ctx);

    h = _fq_vec_init(2 * n - 1, ctx);
    t = _fq_vec_init(2 * n - 1, ctx);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _fq_vec_set(B->rows[i], poly1 + i * m, m, ctx);

    _fq_vec_set(B->rows[i], poly1 + i * m, len1 % m, ctx);

    /* Set rows of A to powers of poly2 */
    fq_one(A->rows[0], ctx);
    _fq_vec_set(A->rows[1], poly2, n, ctx);
    tmp = _fq_vec_init(2 * n - 1, ctx);
    for (i = 2; i < m; i++)
    {
        _fq_poly_mulmod_preinv(tmp, A->rows[i - 1], n, poly2, n, poly3, len3,
                               poly3inv, len3inv, ctx);
        _fq_vec_set(A->rows[i], tmp, n, ctx);
    }
    _fq_vec_clear(tmp, 2 * n - 1, ctx);

    fq_mat_mul(C, B, A, ctx);

    /* Evaluate block composition using the Horner scheme */
    _fq_vec_set(res, C->rows[m - 1], n, ctx);
    _fq_poly_mulmod_preinv(h, A->rows[m - 1], n, poly2, n, poly3, len3,
                           poly3inv, len3inv, ctx);

    for (i = m - 2; i >= 0; i--)
    {
        _fq_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                               poly3inv, len3inv, ctx);
        _fq_poly_add(res, t, n, C->rows[i], n, ctx);
    }

    _fq_vec_clear(h, 2 * n - 1, ctx);
    _fq_vec_clear(t, 2 * n - 1, ctx);

    fq_mat_clear(A, ctx);
    fq_mat_clear(B, ctx);
    fq_mat_clear(C, ctx);
}

void
fq_poly_compose_mod_brent_kung_preinv(fq_poly_t res, const fq_poly_t poly1,
                                      const fq_poly_t poly2,
                                      const fq_poly_t poly3,
                                      const fq_poly_t poly3inv,
                                      const fq_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len3inv = poly3inv->length;
    slong len = len3 - 1;
    slong vec_len = FLINT_MAX(len3 - 1, len2);

    fq_struct *ptr2;
    fq_t inv3;

    if (len3 == 0)
    {
        printf("Exception: division by zero in "
               "fq_poly_compose_mod_brent_kung_preinv\n");
        abort();
    }

    if (len1 >= len3)
    {
        printf("Exception: fq_poly_compose_brent_kung: the degree of the"
               " first polynomial must be smaller than that of the modulus\n");
        abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        fq_poly_zero(res, ctx);
        return;
    }

    if (len1 == 1)
    {
        fq_poly_set(res, poly1, ctx);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        fq_poly_t tmp;
        fq_poly_init(tmp, ctx);
        fq_poly_compose_mod_brent_kung_preinv(tmp, poly1, poly2, poly3,
                                              poly3inv, ctx);
        fq_poly_swap(tmp, res, ctx);
        fq_poly_clear(tmp, ctx);
        return;
    }

    ptr2 = _fq_vec_init(vec_len, ctx);

    if (len2 <= len)
    {
        _fq_vec_set(ptr2, poly2->coeffs, len2, ctx);
        _fq_vec_zero(ptr2 + len2, vec_len - len2, ctx);
    }
    else
    {
        fq_init(inv3, ctx);
        fq_inv(inv3, poly3->coeffs + len, ctx);
        _fq_poly_rem(ptr2, poly2->coeffs, len2,
                     poly3->coeffs, len3, inv3, ctx);
        fq_clear(inv3, ctx);
    }

    fq_poly_fit_length(res, len, ctx);
    _fq_poly_compose_mod_brent_kung_preinv(res->coeffs, poly1->coeffs, len1,
                                           ptr2, poly3->coeffs, len3,
                                           poly3inv->coeffs, len3inv, ctx);
    _fq_poly_set_length(res, len, ctx);
    _fq_poly_normalise(res, ctx);

    _fq_vec_clear(ptr2, vec_len, ctx);
}
