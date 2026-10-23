/*
    Copyright (C) 2011, 2025 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
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
_gr_poly_preinv_compose_mod_brent_kung(
    gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2,
    const gr_poly_preinv_t P,
    gr_ctx_t ctx)
{
    gr_mat_t A, B, C;
    gr_ptr h;
    slong i, n, m;
    slong len3 = P->lenf;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    n = len3 - 1;

    if (len3 == 1)
        return status;

    if (len1 == 1)
        return gr_set(res, poly1, ctx);

    if (len3 == 2)
        return _gr_poly_evaluate(res, poly1, len1, poly2, ctx);

    /* Limitation of Brent-Kung */
    if (len1 >= len3)
        return GR_UNABLE;

    m = n_sqrt(n) + 1;

    gr_mat_init(A, m, n, ctx);
    gr_mat_init(B, m, m, ctx);
    gr_mat_init(C, m, n, ctx);

    GR_TMP_INIT_VEC(h, n, ctx);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        status |= _gr_vec_set(gr_mat_entry_ptr(B, i, 0, ctx), GR_ENTRY(poly1, i * m, sz), m, ctx);

    status |= _gr_vec_set(gr_mat_entry_ptr(B, i, 0, ctx), GR_ENTRY(poly1, i * m, sz), len1 % m, ctx);

    /* Set rows of A to powers of poly2 */
    status |= gr_one(gr_mat_entry_ptr(A, 0, 0, ctx), ctx);
    status |= _gr_vec_set(gr_mat_entry_ptr(A, 1, 0, ctx), poly2, n, ctx);

    for (i = 2; i < m; i++)
    {
        /* Assume that squaring is better. XXX: this depends on the ring. */
#if 1
        status |= _gr_poly_preinv_mulmod(gr_mat_entry_ptr(A, i, 0, ctx),
                    gr_mat_entry_srcptr(A, (i + 1) / 2, 0, ctx), n,
                    gr_mat_entry_srcptr(A, i / 2, 0, ctx), n, P, ctx);
#else
        status |= _gr_poly_preinv_mulmod(gr_mat_entry_ptr(A, i, 0, ctx),
                    gr_mat_entry_srcptr(A, i - 1, 0, ctx), n,
                    poly2, n, P, ctx);
#endif
    }

    status |= gr_mat_mul(C, B, A, ctx);

    /* Evaluate block composition */
    status |= _gr_poly_preinv_mulmod(h, gr_mat_entry_srcptr(A, m - 1, 0, ctx), n, poly2, n, P, ctx);
    status |= _gr_poly_preinv_mod_matrix_rows_evaluate(res, C, h, n, P, ctx);

    GR_TMP_CLEAR_VEC(h, n, ctx);

    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);

    return status;
}

int
_gr_poly_compose_mod_brent_kung_preinv(
    gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2,
    gr_srcptr poly3, slong len3,
    gr_srcptr poly3inv, slong len3inv,
    gr_ctx_t ctx)
{
    gr_poly_preinv_t P;
    _gr_poly_preinv_init_newton_shallow(P, poly3, len3, poly3inv, len3inv, ctx);
    return _gr_poly_preinv_compose_mod_brent_kung(res, poly1, len1, poly2, P, ctx);
}

int
gr_poly_compose_mod_brent_kung_preinv(gr_poly_t res,
                                      const gr_poly_t poly1,
                                      const gr_poly_t poly2,
                                      const gr_poly_t poly3,
                                      const gr_poly_t poly3inv,
                                      gr_ctx_t ctx)
{
    return gr_poly_compose_mod_preinv_wrapper(_gr_poly_compose_mod_brent_kung_preinv, res, poly1, poly2, poly3, poly3inv, ctx);
}
