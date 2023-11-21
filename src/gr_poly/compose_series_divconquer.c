/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

#include "nmod_poly.h"

int
_gr_poly_compose_series_divconquer(gr_ptr res, gr_srcptr poly1, slong len1,
                            gr_srcptr poly2, slong len2, slong N, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, j, k, n;
    slong *hlen, alloc, powlen;
    gr_ptr v, *h, pow, temp;

    if (len1 <= 2 || len2 == 1)
       return _gr_poly_compose_series_horner(res, poly1, len1, poly2, len2, N, ctx);

    /* Initialisation */

    hlen = (slong *) flint_malloc(((len1 + 1) / 2) * sizeof(slong));

    for (k = 1; (2 << k) < len1; k++) ;

    hlen[0] = hlen[1] = FLINT_MIN(N, ((1 << k) - 1) * (len2 - 1) + 1);
    for (i = k - 1; i > 0; i--)
    {
        slong hi = (len1 + (1 << i) - 1) / (1 << i);
        slong t  = FLINT_MIN(N, ((1 << i) - 1) * (len2 - 1) + 1);
        for (n = (hi + 1) / 2; n < hi; n++)
            hlen[n] = t;
    }
    powlen = FLINT_MIN(N, (1 << k) * (len2 - 1) + 1);

    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += hlen[i];

    GR_TMP_INIT_VEC(v, alloc +  2 * powlen, ctx);
    h = (gr_ptr *) flint_malloc(((len1 + 1) / 2) * sizeof(mp_ptr));
    h[0] = v;
    for (i = 0; i < (len1 - 1) / 2; i++)
    {
        h[i + 1] = GR_ENTRY(h[i], hlen[i], sz);
        hlen[i]  = 0;
    }
    hlen[(len1 - 1) / 2] = 0;

    pow  = GR_ENTRY(v, alloc, sz);
    temp = GR_ENTRY(pow, powlen, sz);

    /* Let's start the actual work */

    /* Note: compose_divconquer checks for zero coefficients here.
       This does not work here if poly1 has been truncated. */
    for (i = 0, j = 0; i < len1 / 2; i++, j += 2)
    {
        status |= _gr_vec_mul_scalar(h[i], poly2, len2, GR_ENTRY(poly1, j + 1, sz), ctx);
        status |= gr_add(h[i], h[i], GR_ENTRY(poly1, j, sz), ctx);
        hlen[i] = len2;
    }

    if ((len1 & WORD(1)))
    {
        status |= gr_set(h[i], GR_ENTRY(poly1, j, sz), ctx);
        hlen[i] = 1;
    }

    powlen = FLINT_MIN(N, 2 * len2 - 1);
    status |= _gr_poly_mullow(pow, poly2, len2, poly2, len2, powlen, ctx);

    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            slong templen = FLINT_MIN(N, powlen + hlen[1] - 1);
            status |= _gr_poly_mullow(temp, pow, powlen, h[1], hlen[1], templen, ctx);
            status |= _gr_poly_add(h[0], temp, templen, h[0], hlen[0], ctx);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }

        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2*i + 1] > 0)
            {
                hlen[i] = FLINT_MIN(N, hlen[2*i + 1] + powlen - 1);
                status |= _gr_poly_mullow(h[i], pow, powlen, h[2*i + 1], hlen[2*i + 1], hlen[i], ctx);
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
            hlen[i] = FLINT_MIN(N, hlen[2*i]);
            status |= _gr_vec_set(h[i], h[2*i], hlen[i], ctx);
        }

        status |= _gr_poly_mullow(temp, pow, powlen, pow, powlen, FLINT_MIN(N, 2 * powlen - 1), ctx);
        powlen = FLINT_MIN(N, 2 * powlen - 1);

        {
            gr_ptr t = pow;
            pow      = temp;
            temp     = t;
        }
    }

    status |= _gr_poly_mullow(res, pow, powlen, h[1], hlen[1],  FLINT_MIN(N, powlen + hlen[1] - 1), ctx);
    status |= _gr_vec_add(res, res, h[0], hlen[0], ctx);

    GR_TMP_CLEAR_VEC(v, alloc +  2 * powlen, ctx);

    flint_free(h);
    flint_free(hlen);

    return status;
}

int
gr_poly_compose_series_divconquer(gr_poly_t res,
                    const gr_poly_t poly1,
                    const gr_poly_t poly2, slong n, gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;
    int status;

    if (len2 != 0)
    {
        truth_t is_zero = gr_is_zero(poly2->coeffs, ctx);

        if (is_zero == T_FALSE)
            return GR_DOMAIN;
        if (is_zero == T_UNKNOWN)
            return GR_UNABLE;
    }

    if (len1 == 0 || n == 0)
        return gr_poly_zero(res, ctx);

    if (len2 == 0 || len1 == 1)
        return gr_poly_set_scalar(res, poly1->coeffs, ctx);

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        gr_poly_fit_length(res, lenr, ctx);
        status = _gr_poly_compose_series_divconquer(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
        _gr_poly_set_length(res, lenr, ctx);
        _gr_poly_normalise(res, ctx);
    }
    else
    {
        gr_poly_t t;
        gr_poly_init2(t, lenr, ctx);
        status = _gr_poly_compose_series_divconquer(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
        _gr_poly_set_length(t, lenr, ctx);
        _gr_poly_normalise(t, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }

    return status;
}
