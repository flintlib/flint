/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* https://pages.cs.wisc.edu/~amos/412/lecture-notes/lecture08.pdf */


/* XXX: clarify n - degree vs length
   switch argument order -> n before basis? */

int
_gr_poly_newton_basis_evaluate(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t, u, v;
    slong i;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len <= 1)
    {
        if (len == 0)
            return gr_zero(res, ctx);
        else
            return gr_set(res, poly, ctx);
    }

    GR_TMP_INIT3(t, u, v, ctx);

    status |= gr_sub(u, x, GR_ENTRY(basis, len - 2, sz), ctx);
    status |= gr_mul(t, u, GR_ENTRY(poly, len - 1, sz), ctx);
    status |= gr_add(t, t, GR_ENTRY(poly, len - 2, sz), ctx);

    for (i = len - 3; i >= 0; i--)
    {
        status |= gr_sub(u, x, GR_ENTRY(basis, i, sz), ctx);
        status |= gr_mul(v, u, t, ctx);
        status |= gr_add(t, v, GR_ENTRY(poly, i, sz), ctx);
    }

    gr_swap(res, t, ctx);

    GR_TMP_CLEAR3(t, u, v, ctx);

    return status;
}

int
gr_poly_newton_basis_evaluate(gr_ptr res, const gr_vec_t basis, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)
{
    if (basis->length < poly->length - 1)
        return GR_UNABLE;

    return _gr_poly_newton_basis_evaluate(res, basis->entries, poly->coeffs, poly->length, x, ctx);
}

int
_gr_poly_newton_basis_interpolate_exact(gr_ptr res, gr_srcptr basis, gr_srcptr ys, slong len, gr_ctx_t ctx)
{
    gr_ptr p, q, t;
    slong i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (res != ys)
        status |= _gr_vec_set(res, ys, len, ctx);

    GR_TMP_INIT3(p, q, t, ctx);

    for (i = 1; i < len; i++)
    {
        status |= gr_set(t, GR_ENTRY(res, i - 1, sz), ctx);

        for (j = i; j < len; j++)
        {
            status |= gr_sub(p, GR_ENTRY(res, j, sz), t, ctx);
            status |= gr_sub(q, GR_ENTRY(basis, j, sz), GR_ENTRY(basis, j - i, sz), ctx);
            gr_swap(t, GR_ENTRY(res, j, sz), ctx);
            status |= gr_divexact(GR_ENTRY(res, j, sz), p, q, ctx);
        }
    }

    GR_TMP_CLEAR3(p, q, t, ctx);

    return status;
}

int
gr_poly_newton_basis_interpolate_exact(gr_poly_t res, const gr_vec_t basis, const gr_vec_t ys, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (basis->length < ys->length)
        return GR_UNABLE;

    gr_poly_fit_length(res, ys->length, ctx);
    status |= _gr_poly_newton_basis_interpolate_exact(res->coeffs, basis->entries, ys->entries, ys->length, ctx);
    _gr_poly_set_length(res, ys->length, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
_gr_poly_newton_basis_interpolate(gr_ptr res, gr_srcptr basis, gr_srcptr ys, slong len, gr_ctx_t ctx)
{
    gr_ptr p, q, t;
    slong i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (res != ys)
        status |= _gr_vec_set(res, ys, len, ctx);

    GR_TMP_INIT3(p, q, t, ctx);

    for (i = 1; i < len; i++)
    {
        status |= gr_set(t, GR_ENTRY(res, i - 1, sz), ctx);

        for (j = i; j < len; j++)
        {
            status |= gr_sub(p, GR_ENTRY(res, j, sz), t, ctx);
            status |= gr_sub(q, GR_ENTRY(basis, j, sz), GR_ENTRY(basis, j - i, sz), ctx);
            gr_swap(t, GR_ENTRY(res, j, sz), ctx);
            status |= gr_div(GR_ENTRY(res, j, sz), p, q, ctx);
        }
    }

    GR_TMP_CLEAR3(p, q, t, ctx);

    return status;
}

int
gr_poly_newton_basis_interpolate(gr_poly_t res, const gr_vec_t basis, const gr_vec_t ys, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (basis->length < ys->length)
        return GR_UNABLE;

    gr_poly_fit_length(res, ys->length, ctx);
    status |= _gr_poly_newton_basis_interpolate(res->coeffs, basis->entries, ys->entries, ys->length, ctx);
    _gr_poly_set_length(res, ys->length, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
_gr_poly_newton_basis_to_monomial(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    slong i, j, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (res != poly)
        status |= _gr_vec_set(res, poly, len, ctx);

    for (i = len - 2; i >= 0; i--)
        for (j = i; j < len - 1; j++)
            status |= gr_submul(GR_ENTRY(poly, j, sz), GR_ENTRY(poly, j + 1, sz), GR_ENTRY(basis, i, sz), ctx);

    return status;
}

int
gr_poly_newton_basis_to_monomial(gr_poly_t res, const gr_vec_t basis, const gr_poly_t poly, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (basis->length < poly->length - 1)
        return GR_UNABLE;

    if (res != poly)
        status |= gr_poly_set(res, poly, ctx);

    status |= _gr_poly_newton_basis_to_monomial(res->coeffs, basis->entries, res->coeffs, poly->length, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
_gr_poly_newton_basis_from_monomial(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    slong i, j, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (res != poly)
        status |= _gr_vec_set(res, poly, len, ctx);

    for (i = 0; i < len - 1; i++)
        for (j = len - 2; j >= i; j--)
            status |= gr_addmul(GR_ENTRY(poly, j, sz), GR_ENTRY(poly, j + 1, sz), GR_ENTRY(basis, i, sz), ctx);

    return status;
}

int
gr_poly_newton_basis_from_monomial(gr_poly_t res, const gr_vec_t basis, const gr_poly_t poly, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (basis->length < poly->length - 1)
        return GR_UNABLE;

    if (res != poly)
        status |= gr_poly_set(res, poly, ctx);

    status |= _gr_poly_newton_basis_from_monomial(res->coeffs, basis->entries, res->coeffs, poly->length, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

