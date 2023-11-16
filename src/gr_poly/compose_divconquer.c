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

/* todo: use TMP temporaries */
int
_gr_poly_compose_divconquer(gr_ptr res, gr_srcptr poly1, slong len1,
                                          gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    slong i, j, k, n;
    slong *hlen, alloc, powlen;
    gr_ptr v, pow, temp;
    gr_ptr * h;
    slong sz;
    int status;

    if (len1 == 1)
        return gr_set(res, poly1, ctx);

    if (len2 == 1)
        return _gr_poly_evaluate(res, poly1, len1, poly2, ctx);

    if (len1 == 2)
        return _gr_poly_compose_horner(res, poly1, len1, poly2, len2, ctx);

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    /* Initialisation */

    hlen = (slong *) flint_malloc(((len1 + 1) / 2) * sizeof(slong));

    for (k = 1; (2 << k) < len1; k++) ;

    hlen[0] = hlen[1] = ((1 << k) - 1) *(len2 - 1) + 1;
    for (i = k - 1; i > 0; i--)
    {
        slong hi = (len1 + (1 << i) - 1) / (1 << i);
        for (n = (hi + 1) / 2; n < hi; n++)
            hlen[n] = ((1 << i) - 1) * (len2 - 1) + 1;
    }
    powlen = (1 << k) * (len2 - 1) + 1;

    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += hlen[i];

    GR_TMP_INIT_VEC(v, alloc + 2 * powlen, ctx);

    h = (gr_ptr *) flint_malloc(((len1 + 1) / 2) * sizeof(gr_ptr));
    h[0] = v;
    for (i = 0; i < (len1 - 1) / 2; i++)
    {
        h[i + 1] = GR_ENTRY(h[i], hlen[i], sz);
        hlen[i] = 0;
    }
    hlen[(len1 - 1) / 2] = 0;

    pow = GR_ENTRY(v, alloc, sz);
    temp = GR_ENTRY(pow, powlen, sz);

    /* Let's start the actual work */

    for (i = 0, j = 0; i < len1 / 2; i++, j += 2)
    {
        if (gr_is_zero(GR_ENTRY(poly1, j + 1, sz), ctx) != T_TRUE)
        {
            status |= _gr_vec_mul_scalar(h[i], poly2, len2, GR_ENTRY(poly1, j + 1, sz), ctx);
            status |= gr_add(h[i], h[i], GR_ENTRY(poly1, j, sz), ctx);
            hlen[i] = len2;
        }
        else if (gr_is_zero(GR_ENTRY(poly1, j, sz), ctx) != T_TRUE)
        {
            status |= gr_set(h[i], GR_ENTRY(poly1, j, sz), ctx);
            hlen[i] = 1;
        }
    }

    if ((len1 & WORD(1)))
    {
        if (gr_is_zero(GR_ENTRY(poly1, j, sz), ctx) != T_TRUE)
        {
            status |= gr_set(h[i], GR_ENTRY(poly1, j, sz), ctx);
            hlen[i] = 1;
        }
    }

    status |= _gr_poly_mul(pow, poly2, len2, poly2, len2, ctx);
    powlen = 2 * len2 - 1;

    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            slong templen = powlen + hlen[1] - 1;
            status |= _gr_poly_mul(temp, pow, powlen, h[1], hlen[1], ctx);
            status |= _gr_poly_add(h[0], temp, templen, h[0], hlen[0], ctx);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }

        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2*i + 1] > 0)
            {
                status |= _gr_poly_mul(h[i], pow, powlen, h[2*i + 1], hlen[2*i + 1], ctx);
                hlen[i] = hlen[2*i + 1] + powlen - 1;
            }
            else
            {
                hlen[i] = 0;
            }

            status |= _gr_poly_add(h[i], h[i], hlen[i], h[2*i], hlen[2*i], ctx);
            hlen[i] = FLINT_MAX(hlen[i], hlen[2*i]);
        }
        if ((n & WORD(1)))
        {
            status |= _gr_vec_set(h[i], h[2*i], hlen[2*i], ctx);
            hlen[i] = hlen[2*i];
        }

        status |= _gr_poly_mul(temp, pow, powlen, pow, powlen, ctx);
        powlen += powlen - 1;

        {
            gr_ptr t = pow;
            pow = temp;
            temp = t;
        }
    }

    status |= _gr_poly_mul(res, pow, powlen, h[1], hlen[1], ctx);
    status |= _gr_vec_add(res, res, h[0], hlen[0], ctx);

    GR_TMP_CLEAR_VEC(v, alloc + 2 * powlen, ctx);

    flint_free(h);
    flint_free(hlen);

    return status;
}

int gr_poly_compose_divconquer(gr_poly_t res,
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
            status = _gr_poly_compose_divconquer(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
        }
        else
        {
            gr_poly_t t;
            gr_poly_init2(t, lenr, ctx);
            status = _gr_poly_compose_divconquer(t->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
            gr_poly_swap(res, t, ctx);
            gr_poly_clear(t, ctx);
        }

        _gr_poly_set_length(res, lenr, ctx);
        _gr_poly_normalise(res, ctx);
        return status;
    }
}
