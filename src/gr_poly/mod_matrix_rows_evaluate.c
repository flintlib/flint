/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2025, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"

/*
    Given a matrix A with r rows, each of length n = len3 - 1 representing
    a polynomial c_i, computes res = sum_i c_i h^i mod poly3, where
    h is a polynomial of length n (with reduction using the preinverted
    reciprocal of poly3). This is the final step of Brent-Kung modular
    composition.
*/

static int
_gr_poly_preinv_mod_matrix_rows_evaluate_horner(gr_ptr res, const gr_mat_t A,
    gr_srcptr h, slong n, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    gr_ptr t;
    slong len = A->r;
    slong i;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(t, n, ctx);

    status |= _gr_vec_set(res, gr_mat_entry_srcptr(A, len - 1, 0, ctx), n, ctx);

    for (i = len - 2; i >= 0; i--)
    {
        status |= _gr_poly_preinv_mulmod(t, res, n, h, n, P, ctx);
        status |= _gr_poly_add(res, t, n, gr_mat_entry_srcptr(A, i, 0, ctx), n, ctx);
    }

    GR_TMP_CLEAR_VEC(t, n, ctx);

    return status;
}

/*
    Rectangular splitting: with m ~ sqrt(len), precompute h^0, ..., h^m
    and evaluate blocks of m coefficients using unreduced products,
    performing only one reduction per block. This reduces the number of
    modular reductions from len to about len/m + m.
*/
static int
_gr_poly_preinv_mod_matrix_rows_evaluate_rectangular(gr_ptr res, const gr_mat_t A,
    gr_srcptr h, slong n, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    gr_ptr xs, s, t, q, u;
    slong len = A->r;
    slong i, j, m, r;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    m = n_sqrt(len) + 1;
    m = FLINT_MIN(m, 30);
    r = (len + m - 1) / m;

    GR_TMP_INIT_VEC(xs, (m + 1) * n, ctx);
    GR_TMP_INIT_VEC(s, 8 * n, ctx);
    t = GR_ENTRY(s, 2 * n, sz);
    q = GR_ENTRY(s, 4 * n, sz);
    u = GR_ENTRY(s, 6 * n, sz);

#define XP(ii) GR_ENTRY(xs, (ii) * n, sz)
#define COEFF(ii) gr_mat_entry_srcptr(A, (ii), 0, ctx)

    /* Compute powers of h */
    for (i = 0; i <= m; i++)
    {
        if (i == 0)
        {
            status |= gr_one(XP(0), ctx);
            status |= _gr_vec_zero(GR_ENTRY(XP(0), 1, sz), n - 1, ctx);
        }
        else if (i == 1)
        {
            status |= _gr_vec_set(XP(1), h, n, ctx);
        }
        else
        {
            status |= _gr_poly_preinv_mulmod(XP(i), XP(i / 2), n,
                XP((i + 1) / 2), n, P, ctx);
        }
    }

    status |= _gr_vec_set(s, COEFF((r - 1) * m), n, ctx);
    status |= _gr_vec_zero(GR_ENTRY(s, n, sz), n - 1, ctx);
    for (j = 0; j < len - (r - 1) * m - 1; j++)
    {
        status |= _gr_poly_mul(t, XP(1 + j), n, COEFF((r - 1) * m + 1 + j), n, ctx);
        status |= _gr_vec_add(s, s, t, 2 * n - 1, ctx);
    }
    status |= _gr_poly_preinv_divrem(q, res, s, 2 * n - 1, P, ctx);

    for (i = r - 2; i >= 0; i--)
    {
        status |= _gr_vec_set(s, COEFF(i * m), n, ctx);
        status |= _gr_vec_zero(GR_ENTRY(s, n, sz), n - 1, ctx);
        for (j = 1; j < m; j++)
        {
            status |= _gr_poly_mul(t, XP(j), n, COEFF(i * m + j), n, ctx);
            status |= _gr_vec_add(s, s, t, 2 * n - 1, ctx);
        }
        status |= _gr_poly_preinv_divrem(q, u, s, 2 * n - 1, P, ctx);
        status |= _gr_poly_preinv_mulmod(t, res, n, XP(m), n, P, ctx);
        status |= _gr_vec_add(res, t, u, n, ctx);
    }

#undef XP
#undef COEFF

    GR_TMP_CLEAR_VEC(xs, (m + 1) * n, ctx);
    GR_TMP_CLEAR_VEC(s, 8 * n, ctx);

    return status;
}

int
_gr_poly_preinv_mod_matrix_rows_evaluate(gr_ptr res, const gr_mat_t A,
    gr_srcptr h, slong n, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    FLINT_ASSERT(A->c == n);
    FLINT_ASSERT(n == P->lenf - 1);

    if (A->r <= 10)
        return _gr_poly_preinv_mod_matrix_rows_evaluate_horner(res, A, h, n, P, ctx);
    else
        return _gr_poly_preinv_mod_matrix_rows_evaluate_rectangular(res, A, h, n, P, ctx);
}

int
_gr_poly_mod_matrix_rows_evaluate(gr_ptr res, const gr_mat_t A,
    gr_srcptr h, slong n, gr_srcptr poly3, slong len3,
    gr_srcptr poly3inv, slong len3inv, gr_ctx_t ctx)
{
    gr_poly_preinv_t P;
    _gr_poly_preinv_init_newton_shallow(P, poly3, len3, poly3inv, len3inv, ctx);
    return _gr_poly_preinv_mod_matrix_rows_evaluate(res, A, h, n, P, ctx);
}
