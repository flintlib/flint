/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_mulmod(gr_ptr res, gr_srcptr poly1, slong len1,
                           gr_srcptr poly2, slong len2,
                           gr_srcptr f, slong lenf,
                           gr_ctx_t ctx)
{
    gr_ptr T, Q;
    slong lenT, lenQ;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    lenT = len1 + len2 - 1;
    lenQ = lenT - lenf + 1;

    GR_TMP_INIT_VEC(T, lenT + lenQ, ctx);
    Q = GR_ENTRY(T, lenT, sz);

    status |= _gr_poly_mul(T, poly1, len1, poly2, len2, ctx);
    status |= _gr_poly_divrem(Q, res, T, lenT, f, lenf, ctx);

    GR_TMP_CLEAR_VEC(T, lenT + lenQ, ctx);

    return status;
}

int
gr_poly_mulmod(gr_poly_t res,
                          const gr_poly_t poly1,
                          const gr_poly_t poly2,
                          const gr_poly_t f,
                          gr_ctx_t ctx)
{
    slong len1, len2, lenf;
    gr_ptr fcoeffs;
    int status = GR_SUCCESS;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf == 0)
        return GR_DOMAIN;

    if (lenf == 1 || len1 == 0 || len2 == 0)
        return gr_poly_zero(res, ctx);

    if (len1 + len2 - lenf > 0)
    {
        if (f == res)
        {
            GR_TMP_INIT_VEC(fcoeffs, lenf, ctx);
            status |= _gr_vec_set(fcoeffs, f->coeffs, lenf, ctx);
        }
        else
            fcoeffs = f->coeffs;

        gr_poly_fit_length(res, len1 + len2 - 1, ctx);
        status |= _gr_poly_mulmod(res->coeffs,
                                   poly1->coeffs, len1,
                                   poly2->coeffs, len2, fcoeffs, lenf, ctx);

        if (f == res)
            GR_TMP_CLEAR_VEC(fcoeffs, lenf, ctx);

        _gr_poly_set_length(res, lenf - 1, ctx);
        _gr_poly_normalise(res, ctx);
        return status;
    }
    else
    {
        return gr_poly_mul(res, poly1, poly2, ctx);
    }
}
