/*
    Copyright (C) 2023, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_mulmid_generic(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nhi != 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

    /* Todo: a reasonable generic implementation. Currently many rings
       implement a good mullow, so fall back on that. */
    if (ctx->methods[GR_METHOD_POLY_MULLOW] != (gr_funcptr) _gr_poly_mullow_generic)
    {
        if (nlo == 0)
            return _gr_poly_mullow(res, poly1, len1, poly2, len2, nhi, ctx);

        gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);
        gr_ptr tmp;
        slong i, sz = ctx->sizeof_elem;
        int status;

        /* Only the low part of the temporary product needs new elements:
           the high part borrows the elements of res (moved shallowly back
           and forth) so that their storage is reused. */
        tmp = GR_TMP_ALLOC(nhi * sz);
        _gr_vec_init(tmp, nlo, ctx);
        for (i = 0; i < nhi - nlo; i++)
            set_shallow(GR_ENTRY(tmp, nlo + i, sz), GR_ENTRY(res, i, sz), ctx);

        status = _gr_poly_mullow(tmp, poly1, len1, poly2, len2, nhi, ctx);

        for (i = 0; i < nhi - nlo; i++)
            set_shallow(GR_ENTRY(res, i, sz), GR_ENTRY(tmp, nlo + i, sz), ctx);
        _gr_vec_clear(tmp, nlo, ctx);
        GR_TMP_FREE(tmp, nhi * sz);

        return status;
    }
    else
    {
        return _gr_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, ctx);
    }
}

int
gr_poly_mulmid(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong nlo, slong nhi, gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    int status;
    slong len;

    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi >= 0);

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
        return gr_poly_zero(res, ctx);

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, len, ctx);
        status = _gr_poly_mulmid(t->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len, ctx);
        status = _gr_poly_mulmid(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi, ctx);
    }

    _gr_poly_set_length_normalise(res, len, ctx);
    return status;
}

