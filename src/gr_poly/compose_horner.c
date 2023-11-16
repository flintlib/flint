/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_compose_horner(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len1 == 1)
    {
        return gr_set(res, poly1, ctx);
    }
    else if (len2 == 1)
    {
        return _gr_poly_evaluate(res, poly1, len1, poly2, ctx);
    }
    else if (len1 == 2)
    {
        slong sz = ctx->sizeof_elem;
        status |= _gr_vec_mul_scalar(res, poly2, len2, GR_ENTRY(poly1, 1, sz), ctx);
        status |= gr_add(res, res, poly1, ctx);
        return status;
    }
    else
    {
        slong alloc = (len1 - 1) * (len2 - 1) + 1;
        slong i = len1 - 1, lenr = len2;
        slong sz = ctx->sizeof_elem;
        gr_ptr t, t1, t2;

        GR_TMP_INIT_VEC(t, alloc, ctx);

        if (len1 % 2 == 0)
        {
            t1 = res;
            t2 = t;
        }
        else
        {
            t1 = t;
            t2 = res;
        }

        /* Perform the first two steps as one,
            "res = a(m) * poly2 + a(m-1)". */
        status |= _gr_vec_mul_scalar(t1, poly2, len2, GR_ENTRY(poly1, i, sz), ctx);
        i--;
        status |= gr_add(t1, t1, GR_ENTRY(poly1, i, sz), ctx);

        while (i--)
        {
            status |= _gr_poly_mul(t2, t1, lenr, poly2, len2, ctx);
            lenr += len2 - 1;

            {
                void *t_ = t1;
                t1 = t2;
                t2 = t_;
            }

            status |= gr_add(t1, t1, GR_ENTRY(poly1, i, sz), ctx);
        }

        GR_TMP_CLEAR_VEC(t, alloc, ctx);

        return status;
    }
}

int gr_poly_compose_horner(gr_poly_t res,
    const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0)
    {
        return gr_poly_zero(res, ctx);
    }
    else if (len1 == 1 || len2 == 0)
    {
        return gr_poly_set_scalar(res, poly1->coeffs, ctx);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;
        int status;

        if (res != poly1 && res != poly2)
        {
            gr_poly_fit_length(res, lenr, ctx);
            status = _gr_poly_compose_horner(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
        }
        else
        {
            gr_poly_t t;
            gr_poly_init2(t, lenr, ctx);
            status = _gr_poly_compose_horner(t->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
            gr_poly_swap(res, t, ctx);
            gr_poly_clear(t, ctx);
        }

        _gr_poly_set_length(res, lenr, ctx);
        _gr_poly_normalise(res, ctx);
        return status;
    }
}
