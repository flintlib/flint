/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "thread_support.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "factor_impl.h"

typedef struct
{
    gr_poly_struct * res;
    const gr_mat_struct * C;
    gr_srcptr h;
    slong n, k;
    const gr_poly_preinv_struct * P;
    gr_ctx_struct * ctx;
    int * status;
}
vec_eval_args_t;

static void
_vec_eval_worker(slong j, void * arg)
{
    vec_eval_args_t * args = arg;
    gr_mat_t Cw;
    gr_mat_window_init(Cw, args->C, j * args->k, 0, (j + 1) * args->k, args->n, args->ctx);
    args->status[j] = _gr_poly_preinv_mod_matrix_rows_evaluate(args->res[j].coeffs, Cw, args->h, args->n, args->P, args->ctx);
    gr_mat_window_clear(Cw, args->ctx);
}

/*
    Sets res[j] = polys[j](g) mod poly for 0 <= j < l, where the polys
    have length less than len = len(poly) and g has length glen <= len - 1
    (padded with zeros). Uses a single matrix product for all polynomials
    with rectangular splitting parameter m ~ sqrt(n l), which reduces the
    number of modular multiplications compared to l separate compositions.
    Each res[j] must have space for len - 1 coefficients.
*/
int
_gr_poly_preinv_compose_mod_brent_kung_vec(gr_poly_struct * res,
    const gr_poly_struct * polys, slong FLINT_UNUSED(lenpolys), slong l,
    gr_srcptr g, slong glen, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong len = P->lenf;
    gr_mat_t A, B, C;
    gr_ptr h;
    slong i, j, k, n, m, len1;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    n = len - 1;
    m = n_sqrt(n * l) + 1;
    k = len / m + 1;

    gr_mat_init(A, m, n, ctx);
    gr_mat_init(B, k * l, m, ctx);
    gr_mat_init(C, k * l, n, ctx);
    GR_TMP_INIT_VEC(h, n, ctx);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < l; j++)
    {
        len1 = (polys + j)->length;

        for (i = 0; i < len1 / m; i++)
            status |= _gr_vec_set(gr_mat_entry_ptr(B, i + j * k, 0, ctx), GR_ENTRY((polys + j)->coeffs, i * m, sz), m, ctx);

        status |= _gr_vec_set(gr_mat_entry_ptr(B, i + j * k, 0, ctx), GR_ENTRY((polys + j)->coeffs, i * m, sz), len1 % m, ctx);
    }

    /* Set rows of A to powers of g */
    status |= gr_one(gr_mat_entry_ptr(A, 0, 0, ctx), ctx);
    if (m > 1)
    {
        status |= _gr_vec_set(gr_mat_entry_ptr(A, 1, 0, ctx), g, glen, ctx);
        status |= _gr_vec_zero(gr_mat_entry_ptr(A, 1, glen, ctx), n - glen, ctx);
    }
    for (i = 2; i < m; i++)
        status |= _gr_poly_preinv_mulmod(gr_mat_entry_ptr(A, i, 0, ctx),
                gr_mat_entry_srcptr(A, (i + 1) / 2, 0, ctx), n,
                gr_mat_entry_srcptr(A, i / 2, 0, ctx), n, P, ctx);

    status |= gr_mat_mul(C, B, A, ctx);

    /* h = g^m */
    status |= _gr_poly_preinv_mulmod(h, gr_mat_entry_srcptr(A, m - 1, 0, ctx), n,
                gr_mat_entry_srcptr(A, 1, 0, ctx), n, P, ctx);

    /* Evaluate the block compositions (independent for each polynomial:
       in parallel when the ring is threadsafe) */
    if (status == GR_SUCCESS)
    {
        vec_eval_args_t args;
        args.res = res;
        args.C = C;
        args.h = h;
        args.n = n;
        args.k = k;
        args.P = P;
        args.ctx = ctx;
        args.status = flint_calloc(l, sizeof(int));

        if (l > 1 && n > gr_poly_factor_threaded_cutoff && flint_get_num_available_threads() > 1 && gr_ctx_is_threadsafe(ctx) == T_TRUE)
            flint_parallel_do(_vec_eval_worker, &args, l, -1, FLINT_PARALLEL_UNIFORM);
        else
            for (j = 0; j < l; j++)
                _vec_eval_worker(j, &args);

        for (j = 0; j < l; j++)
            status |= args.status[j];

        flint_free(args.status);
    }

    GR_TMP_CLEAR_VEC(h, n, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);

    return status;
}

int
_gr_poly_compose_mod_brent_kung_vec_preinv(gr_poly_struct * res,
    const gr_poly_struct * polys, slong lenpolys, slong l,
    gr_srcptr g, slong glen, gr_srcptr poly, slong len,
    gr_srcptr polyinv, slong leninv, gr_ctx_t ctx)
{
    gr_poly_preinv_t P;
    _gr_poly_preinv_init_newton_shallow(P, poly, len, polyinv, leninv, ctx);
    return _gr_poly_preinv_compose_mod_brent_kung_vec(res, polys, lenpolys, l, g, glen, P, ctx);
}

int
gr_poly_preinv_compose_mod_brent_kung_vec(gr_poly_struct * res,
    const gr_poly_struct * polys, slong len1, slong n,
    const gr_poly_t g, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong len2 = P->lenf;
    slong len3, i;
    int status = GR_SUCCESS;

    if (n > len1)
        return GR_DOMAIN;

    for (i = 0; i < n; i++)
    {
        len3 = (polys + i)->length;
        if (len3 >= len2)
            return GR_UNABLE;
    }

    if (n == 0)
        return GR_SUCCESS;

    if (len2 == 0)
        return GR_DOMAIN;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            status |= gr_poly_zero(res + i, ctx);
        return status;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            status |= gr_poly_set(res + i, polys + i, ctx);
        return status;
    }

    if (g->length > len2 - 1)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status |= gr_poly_preinv_rem(t, g, P, ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_preinv_compose_mod_brent_kung_vec(res, polys, len1, n, t, P, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    for (i = 0; i < n; i++)
    {
        if (res + i == g)
            return GR_UNABLE;  /* aliasing not supported */
        gr_poly_fit_length(res + i, len2 - 1, ctx);
        _gr_poly_set_length(res + i, len2 - 1, ctx);
    }

    status |= _gr_poly_preinv_compose_mod_brent_kung_vec(res, polys, len1, n, g->coeffs, g->length, P, ctx);

    for (i = 0; i < n; i++)
        _gr_poly_normalise(res + i, ctx);

    return status;
}

int
gr_poly_compose_mod_brent_kung_vec_preinv(gr_poly_struct * res,
    const gr_poly_struct * polys, slong len1, slong n,
    const gr_poly_t g, const gr_poly_t poly, const gr_poly_t polyinv, gr_ctx_t ctx)
{
    slong len2 = poly->length;
    slong len3, i;
    int status = GR_SUCCESS;

    if (n > len1)
        return GR_DOMAIN;

    for (i = 0; i < n; i++)
    {
        len3 = (polys + i)->length;

        /* Limitation of Brent-Kung */
        if (len3 >= len2)
            return GR_UNABLE;
    }

    if (n == 0)
        return GR_SUCCESS;

    if (len2 == 0)
        return GR_DOMAIN;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            status |= gr_poly_zero(res + i, ctx);
        return status;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            status |= gr_poly_set(res + i, polys + i, ctx);
        return status;
    }

    if (g->length > len2 - 1)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status |= gr_poly_rem(t, g, poly, ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_compose_mod_brent_kung_vec_preinv(res, polys, len1, n, t, poly, polyinv, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    for (i = 0; i < n; i++)
    {
        if (res + i == g || res + i == poly || res + i == polyinv)
            return GR_UNABLE;  /* aliasing not supported */
        gr_poly_fit_length(res + i, len2 - 1, ctx);
        _gr_poly_set_length(res + i, len2 - 1, ctx);
    }

    status |= _gr_poly_compose_mod_brent_kung_vec_preinv(res, polys, len1, n,
                g->coeffs, g->length, poly->coeffs, len2, polyinv->coeffs,
                polyinv->length, ctx);

    for (i = 0; i < n; i++)
        _gr_poly_normalise(res + i, ctx);

    return status;
}
