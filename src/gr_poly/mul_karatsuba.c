/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* TODO: a self-recursive variant that can reuse scratch space. */
int
_gr_poly_mul_karatsuba(gr_ptr res, gr_srcptr f, slong flen, gr_srcptr g, slong glen, gr_ctx_t ctx)
{
    slong m, f1len, g1len, tlen, ulen, vlen, alloc;
    slong sz = ctx->sizeof_elem;
    gr_ptr t, u, v;
    gr_srcptr f0, f1, g0, g1;
    int squaring = (f == g) && (flen == glen);
    int status = GR_SUCCESS;

    FLINT_ASSERT(flen >= 0);
    FLINT_ASSERT(glen >= 0);
    FLINT_ASSERT(res != f);
    FLINT_ASSERT(res != g);

    /* TODO: should explicitly call basecase mul. */
    if (flen == 1 || glen == 1)
        return _gr_poly_mullow_generic(res, f, flen, g, glen, flen + glen - 1, ctx);

    /* Split at X = x^m */
    /* res = f0 g0 + (f0 g1 + f1 g0) X + f1 g1 X^2
           = f0 g0 + ((f0 + f1) (g0 + g1) - f0 g0 - f1 g1) X + f1 g1 X^2 */
    m = (FLINT_MIN(flen, glen) + 1) / 2;

    f0 = f;
    g0 = g;
    f1 = GR_ENTRY(f, m, sz);
    g1 = GR_ENTRY(g, m, sz);
    f1len = flen - m;
    g1len = glen - m;

    /* Low part: res[0, ..., 2m-2] = f0 g0. */
    status |= _gr_poly_mul(res, f, m, g, m, ctx);
    /* res[2m-1] = 0 (could be avoided if we split the summation later, which is more annoying) */
    status |= gr_zero(GR_ENTRY(res, 2 * m - 1, sz), ctx);

    /* High part: res[2m, ..., flen+glen-2] = f1 g1. */
    status |= _gr_poly_mul(GR_ENTRY(res, 2 * m, sz), f1, f1len, g1, g1len, ctx);

    /* Temporary space for the middle part. */
    tlen = FLINT_MAX(m, f1len);
    ulen = FLINT_MAX(m, g1len);
    vlen = tlen + ulen - 1;
    alloc = tlen + (squaring ? 0 : ulen) + vlen;

    GR_TMP_INIT_VEC(t, alloc, ctx);
    u = GR_ENTRY(t, tlen, sz);
    v = squaring ? u : GR_ENTRY(u, ulen, sz);

    /* t = f0 + f1 */
    status |= _gr_poly_add(t, f0, m, f1, f1len, ctx);

    if (squaring)
    {
        status |= _gr_poly_mul(v, t, tlen, t, tlen, ctx);
    }
    else
    {
        /* u = g0 + g1 */
        status |= _gr_poly_add(u, g0, m, g1, g1len, ctx);
        /* v = (f0 + f1) (g0 + g1) */
        status |= _gr_poly_mul(v, t, tlen, u, ulen, ctx);
    }

    /* v -= f0 g0 */
    status |= _gr_vec_sub(v, v, res, 2 * m - 1, ctx);
    /* v -= f1 g1 */
    status |= _gr_vec_sub(v, v, GR_ENTRY(res, 2 * m, sz), f1len + g1len - 1, ctx);

    /* Finally add the middle part. */
    status |= _gr_poly_add(GR_ENTRY(res, m, sz), GR_ENTRY(res, m, sz), vlen, v, vlen, ctx);

    GR_TMP_CLEAR_VEC(t, alloc, ctx);

    return status;
}

int
gr_poly_mul_karatsuba(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    slong len_out;
    int status;

    if (poly1->length == 0 || poly2->length == 0)
        return gr_poly_zero(res, ctx);

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, len_out, ctx);
        status = _gr_poly_mul_karatsuba(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len_out, ctx);
        status = _gr_poly_mul_karatsuba(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, ctx);
    }

    _gr_poly_set_length(res, len_out, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
