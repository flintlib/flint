/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"

void
_nmod_poly_compose_mod_brent_kung_preinv(mp_ptr res, mp_srcptr poly1,
                  slong len1, mp_srcptr poly2, mp_srcptr poly3, slong len3,
                                 mp_srcptr poly3inv, slong len3inv, nmod_t mod)
{
    nmod_mat_t A, B, C;
    mp_ptr t, h;
    slong i, n, m;

    n = len3 - 1;

    if (len3 == 1)
        return;

    if (len1 == 1)
    {
        res[0] = poly1[0];
        return;
    }

    if (len3 == 2)
    {
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
        return;
    }

    m = n_sqrt(n) + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, n, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1/m; i++)
        _nmod_vec_set(B->rows[i], poly1 + i*m, m);

    _nmod_vec_set(B->rows[i], poly1 + i*m, len1%m);

    /* Set rows of A to powers of poly2 */
    _nmod_poly_powers_mod_preinv_naive(A->rows, poly2, n,
                                       m, poly3, len3, poly3inv, len3inv, mod);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res, C->rows[m - 1], n);
    _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, poly2, n,
                                           poly3, len3, poly3inv, len3inv,mod);

    for (i = m - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                                                       poly3inv, len3inv, mod);
        _nmod_poly_add(res, t, n, C->rows[i], n, mod);
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                           const nmod_poly_t poly3, const nmod_poly_t poly3inv)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    mp_ptr ptr2;

    if (len3 == 0)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_compose_mod_brent_kung_preinv): Division by zero.\n");
    }

    if (len1 >= len3)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_compose_mod_brent_kung_preinv): "
                "The degree of the first polynomial must be smaller than that of the modulus.\n");
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1 || res == poly3inv)
    {
        nmod_poly_t tmp;
        nmod_poly_init_mod(tmp, res->mod);
        nmod_poly_compose_mod_brent_kung_preinv(tmp, poly1, poly2,
                                                              poly3, poly3inv);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    ptr2 = _nmod_vec_init(len);

    if (len2 <= len)
    {
        flint_mpn_copyi(ptr2, poly2->coeffs, len2);
        flint_mpn_zero(ptr2 + len2, len - len2);
    } else
    {
        _nmod_poly_rem(ptr2, poly2->coeffs, len2,
                             poly3->coeffs, len3, res->mod);
    }

    nmod_poly_fit_length(res, len);

    _nmod_poly_compose_mod_brent_kung_preinv(res->coeffs,
                      poly1->coeffs, len1, ptr2, poly3->coeffs, len3,
                                 poly3inv->coeffs, poly3inv->length, res->mod);

    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}
