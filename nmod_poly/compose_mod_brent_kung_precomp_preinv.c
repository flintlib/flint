/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

void
_nmod_poly_reduce_matrix_mod_poly (nmod_mat_t A, const nmod_mat_t B,
                                   const nmod_poly_t f)
{
    mp_ptr tmp1;
    slong n = f->length - 1;
    slong i, m = n_sqrt(n) + 1;

    nmod_mat_init(A, m, n, f->mod.n);

    A->rows[0][0] = UWORD(1);
    tmp1 = _nmod_vec_init(B->c - f->length + 1);
    for (i= 1; i < m; i++)
        _nmod_poly_divrem (tmp1, A->rows[i], B->rows[i], B->c, f->coeffs,
                           f->length, f->mod);

    _nmod_vec_clear (tmp1);
}

void *
_nmod_poly_precompute_matrix_worker (void * arg_ptr)
{
    nmod_poly_matrix_precompute_arg_t arg =
                           *((nmod_poly_matrix_precompute_arg_t *) arg_ptr);
    /* Set rows of A to powers of poly1 */
    slong i, n, m;

    n = arg.poly2.length - 1;

    m = n_sqrt(n) + 1;

    arg.A.rows[0][0] = UWORD(1);
    _nmod_vec_set(arg.A.rows[1], arg.poly1.coeffs, n);
    for (i = 2; i < m; i++)
        _nmod_poly_mulmod_preinv(arg.A.rows[i], arg.A.rows[i - 1], n,
                                 arg.poly1.coeffs, n, arg.poly2.coeffs, n + 1,
                                 arg.poly2inv.coeffs, n + 1, arg.poly2.mod);
    flint_cleanup();
    return NULL;
}

void
_nmod_poly_precompute_matrix (nmod_mat_t A, mp_srcptr poly1, mp_srcptr poly2,
                              slong len2, mp_srcptr poly2inv, slong len2inv,
                              nmod_t mod)
{
    /* Set rows of A to powers of poly1 */
    slong i, n, m;

    n = len2 - 1;

    m = n_sqrt(n) + 1;

    A->rows[0][0] = UWORD(1);
    _nmod_vec_set(A->rows[1], poly1, n);
    for (i = 2; i < m; i++)
        _nmod_poly_mulmod_preinv(A->rows[i], A->rows[i - 1],
            n, poly1, n, poly2, len2, poly2inv, len2inv, mod);
}

void
nmod_poly_precompute_matrix (nmod_mat_t A, const nmod_poly_t poly1,
                            const nmod_poly_t poly2, const nmod_poly_t poly2inv)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len = len2 - 1;
    slong m = n_sqrt(len) + 1;

    mp_ptr ptr1;

    if (len2 == 0)
    {
        flint_printf("Exception (nmod_poly_precompute_matrix). Division by zero.\n");
        flint_abort();
    }

    if (A->r != m || A->c != len)
    {
        flint_printf("Exception (nmod_poly_precompute_matrix). Wrong dimensions.\n");
        flint_abort();
    }

    if (len2 == 1)
    {
        nmod_mat_zero(A);
        return;
    }

    ptr1 = _nmod_vec_init(len);

    if (len1 <= len)
    {
        flint_mpn_copyi(ptr1, poly1->coeffs, len1);
        flint_mpn_zero(ptr1 + len1, len - len1);
    }
    else
    {
        _nmod_poly_rem(ptr1, poly1->coeffs, len1,
                             poly2->coeffs, len2, A->mod);
    }

    _nmod_poly_precompute_matrix (A, ptr1, poly2->coeffs, len2, poly2inv->coeffs,
                   poly2inv->length, A->mod);

    _nmod_vec_clear (ptr1);
}

void *
_nmod_poly_compose_mod_brent_kung_precomp_preinv_worker(void * arg_ptr)
{
    nmod_poly_compose_mod_precomp_preinv_arg_t arg=
                   *((nmod_poly_compose_mod_precomp_preinv_arg_t*) arg_ptr);
    nmod_mat_t B, C;
    mp_ptr t, h;
    slong i, n, m;

    n = arg.poly3.length - 1;

    if (arg.poly3.length == 1)
    {
        flint_cleanup();
        return NULL;
    }
    if (arg.poly1.length == 1)
    {
        arg.res.coeffs[0] = arg.poly1.coeffs[0];
        flint_cleanup();
        return NULL;
    }

    if (arg.poly3.length == 2)
    {
        arg.res.coeffs[0] = _nmod_poly_evaluate_nmod(arg.poly1.coeffs,
                                             arg.poly1.length, arg.A.rows[1][0],
                                             arg.poly3.mod);
        flint_cleanup();
        return NULL;
    }

    m = n_sqrt(n) + 1;

    nmod_mat_init(B, m, m, arg.poly3.mod.n);
    nmod_mat_init(C, m, n, arg.poly3.mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < arg.poly1.length / m; i++)
        _nmod_vec_set(B->rows[i], arg.poly1.coeffs + i * m, m);

    _nmod_vec_set(B->rows[i], arg.poly1.coeffs + i * m, arg.poly1.length % m);

    nmod_mat_mul(C, B, &arg.A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(arg.res.coeffs, C->rows[m - 1], n);
    _nmod_poly_mulmod_preinv(h, arg.A.rows[m - 1], n, arg.A.rows[1], n,
                             arg.poly3.coeffs, arg.poly3.length,
                             arg.poly3inv.coeffs, arg.poly3inv.length,
                             arg.poly3.mod);

    for (i = m - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod_preinv(t, arg.res.coeffs, n, h, n, arg.poly3.coeffs,
                                 arg.poly3.length, arg.poly3inv.coeffs,
                                 arg.poly3.length, arg.poly3.mod);
        _nmod_poly_add(arg.res.coeffs, t, n, C->rows[i], n, arg.poly3.mod);
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(B);
    nmod_mat_clear(C);
    flint_cleanup();
    return NULL;
}

void
_nmod_poly_compose_mod_brent_kung_precomp_preinv
                            (mp_ptr res, mp_srcptr poly1, slong len1,
                            const nmod_mat_t A, mp_srcptr poly3, slong len3,
                            mp_srcptr poly3inv, slong len3inv, nmod_t mod)
{
    nmod_mat_t B, C;
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
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, A->rows[1][0], mod);
        return;
    }

    m = n_sqrt(n) + 1;

    /* TODO check A*/

    nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, n, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _nmod_vec_set(B->rows[i], poly1 + i*m, m);

    _nmod_vec_set(B->rows[i], poly1 + i*m, len1 % m);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res, C->rows[m - 1], n);
    _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n,
                             poly3, len3, poly3inv, len3inv,mod);

    for (i = m - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                                 poly3inv, len3inv, mod);
        _nmod_poly_add(res, t, n, C->rows[i], n, mod);
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_precomp_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_mat_t A,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv)
{
    slong len1 = poly1->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    if (len3 == 0)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_precomp_preinv). Division by zero.\n");
        flint_abort();
    }

    if (len1 >= len3)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_precomp_preinv). The degree of the \n"
               "first polynomial must be smaller than that of the modulus.\n");
        flint_abort();
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
        nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
        nmod_poly_compose_mod_brent_kung_precomp_preinv(tmp, poly1, A,
                                                poly3, poly3inv);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    nmod_poly_fit_length(res, len);
    _nmod_poly_compose_mod_brent_kung_precomp_preinv(res->coeffs,
        poly1->coeffs, len1, A, poly3->coeffs, len3,
        poly3inv->coeffs, poly3inv->length, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

}
