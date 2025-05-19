/*
    Copyright (C) 2011, 2025 Fredrik Johansson
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"

int
_gr_poly_reduce_matrix_mod_poly(gr_mat_t A,
                                           const gr_mat_t B,
                                           const gr_poly_t f,
                                           gr_ctx_t ctx)
{
    slong n = f->length - 1;
    slong i, m = n_sqrt(n) + 1;
    int status = GR_SUCCESS;

    gr_mat_init(A, m, n, ctx);

    /* todo: use a rem_preinv or at least preinv1 method for divisions */
    /* gr_inv(invf, f->coeffs + (f->length - 1), ctx); */

    status |= gr_one(gr_mat_entry_ptr(A, 0, 0, ctx), ctx);
    for (i = 1; i < m; i++)
        status |= _gr_poly_rem(gr_mat_entry_ptr(A, i, 0, ctx),
                                gr_mat_entry_srcptr(B, i, 0, ctx),
                                B->c, f->coeffs, f->length, ctx);

    return status;
}

int
_gr_poly_precompute_matrix(
    gr_mat_t A,
    gr_srcptr poly1,
    gr_srcptr poly2, slong len2,
    gr_srcptr poly2inv, slong len2inv,
    gr_ctx_t ctx)
{
    /* Set rows of A to powers of poly1 */
    slong i, n, m;
    int status = GR_SUCCESS;

    n = len2 - 1;

    /* XXX: why not just read A->r? */
    m = n_sqrt(n) + 1;

    status |= gr_one(gr_mat_entry_ptr(A, 0, 0, ctx), ctx);
    status |= _gr_vec_set(gr_mat_entry_ptr(A, 1, 0, ctx), poly1, n, ctx);
    for (i = 2; i < m; i++)
        status |= _gr_poly_mulmod_preinv(gr_mat_entry_ptr(A, i, 0, ctx),
                gr_mat_entry_srcptr(A, (i + 1) / 2, 0, ctx), n,
                gr_mat_entry_srcptr(A, i / 2, 0, ctx), n, poly2, len2,
                                          poly2inv, len2inv, ctx);

    return status;
}

int
gr_poly_precompute_matrix(gr_mat_t A,
                                     const gr_poly_t poly1,
                                     const gr_poly_t poly2,
                                     const gr_poly_t poly2inv,
                                     gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len = len2 - 1;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong m = n_sqrt(len) + 1;
    gr_ptr ptr1;

    /* Division by zero */
    if (len2 == 0)
        return GR_DOMAIN;

    /* Wrong dimensions */
    if (A->r != m || A->c != len)
        return GR_DOMAIN;

    if (len2 == 1)
        return gr_mat_zero(A, ctx);

    GR_TMP_INIT_VEC(ptr1, len, ctx);

    if (len1 <= len)
    {
        status |= _gr_vec_set(ptr1, poly1->coeffs, len1, ctx);
        status |= _gr_vec_zero(GR_ENTRY(ptr1, len1, sz), len - len1, ctx);
    }
    else
    {
        status |= _gr_poly_rem(ptr1, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
    }

    status |= _gr_poly_precompute_matrix(A, ptr1, poly2->coeffs, len2,
                                          poly2inv->coeffs, poly2inv->length,
                                          ctx);

    GR_TMP_CLEAR_VEC(ptr1, len, ctx);

    return status;
}

int
_gr_poly_compose_mod_brent_kung_precomp_preinv(
    gr_ptr res,
    gr_srcptr poly1, slong len1,
    const gr_mat_t A,
    gr_srcptr poly3, slong len3,
    gr_srcptr poly3inv, slong len3inv,
    gr_ctx_t ctx)
{
    gr_mat_t B, C;
    gr_ptr t, h;
    slong i, n, m;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    n = len3 - 1;

    if (len3 == 1)
        return status;

    if (len1 == 1)
        return gr_set(res, poly1, ctx);

    if (len3 == 2)
        return _gr_poly_evaluate(res, poly1, len1, gr_mat_entry_srcptr(A, 1, 0, ctx), ctx);

    /* Limitation of Brent-Kung */
    if (len1 >= len3)
        return GR_UNABLE;

    m = n_sqrt(n) + 1;

    /* TODO check A */

    gr_mat_init(B, m, m, ctx);
    gr_mat_init(C, m, n, ctx);

    GR_TMP_INIT_VEC(h, 2 * n, ctx);
    t = GR_ENTRY(h, n, sz);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        status |= _gr_vec_set(gr_mat_entry_ptr(B, i, 0, ctx), GR_ENTRY(poly1, i * m, sz), m, ctx);

    status |= _gr_vec_set(gr_mat_entry_ptr(B, i, 0, ctx), GR_ENTRY(poly1, i * m, sz), len1 % m, ctx);

    status |= gr_mat_mul(C, B, A, ctx);

    /* Evaluate block composition using the Horner scheme */
    status |= _gr_vec_set(res, gr_mat_entry_srcptr(C, m - 1, 0, ctx), n, ctx);
    status |= _gr_poly_mulmod_preinv(h, gr_mat_entry_srcptr(A, m - 1, 0, ctx), n,
                                        gr_mat_entry_srcptr(A, 1, 0, ctx), n,
                                      poly3, len3, poly3inv, len3inv, ctx);

    for (i = m - 2; i >= 0; i--)
    {
        status |= _gr_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                                          poly3inv, len3inv, ctx);
        status |= _gr_poly_add(res, t, n, gr_mat_entry_srcptr(C, i, 0, ctx), n, ctx);
    }

    GR_TMP_CLEAR_VEC(h, 2 * n, ctx);

    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);

    return status;
}

int
gr_poly_compose_mod_brent_kung_precomp_preinv(
    gr_poly_t res,
    const gr_poly_t poly1, const gr_mat_t A,
    const gr_poly_t poly3, const gr_poly_t poly3inv,
    gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;
    int status;

    if (len3 == 0)
        return GR_DOMAIN;

    if (len1 == 0 || len3 == 1)
        return gr_poly_zero(res, ctx);

    if (len1 == 1)
        return gr_poly_set(res, poly1, ctx);

    if (res == poly3 || res == poly1 || res == poly3inv)
    {
        gr_poly_t tmp;
        gr_poly_init(tmp, ctx);
        status = gr_poly_compose_mod_brent_kung_precomp_preinv(tmp, poly1, A,
                                                                 poly3,
                                                                 poly3inv,
                                                                 ctx);
        gr_poly_swap(tmp, res, ctx);
        gr_poly_clear(tmp, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status = _gr_poly_compose_mod_brent_kung_precomp_preinv(res->coeffs,
                                                              poly1->coeffs,
                                                              len1, A,
                                                              poly3->coeffs,
                                                              len3,
                                                              poly3inv->coeffs,
                                                              poly3inv->length,
                                                              ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
