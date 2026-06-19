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

/*
with sign/scale changes to compute the following functions:

    func = 0 -> tan
    func = 1 -> tanh
    func = 2 -> cot
    func = 3 -> coth
    func = 4 -> tan_pi
    func = 5 -> tanh_pi (not yet implemented/used; scalar function missing)
    func = 6 -> cot_pi
    func = 7 -> coth_pi (not yet implemented/used; scalar function missing)
*/
int
_gr_poly_tan_series_sine_cosine(gr_ptr f, gr_srcptr h, slong hlen, slong n, int func, gr_ctx_t ctx)
{
    gr_ptr s, c, t;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    GR_TMP_INIT_VEC(s, 2 * n + 1, ctx);
    c = GR_ENTRY(s, n, sz);
    t = GR_ENTRY(c, n, sz);

    hlen = FLINT_MIN(hlen, n);

    switch (func)
    {
        /* tan(h) = sin(h)/cos(h) */
        case 0:
            status |= _gr_poly_sin_cos_series(s, c, h, hlen, n, ctx);
            status |= _gr_poly_div_series(f, s, n, c, n, n, ctx);
            break;
        /* tanh(h) = -i sin(ih)/cos(ih) */
        case 1:
            status |= gr_i(t, ctx);
            status |= _gr_vec_mul_scalar(s, h, hlen, t, ctx);
            status |= _gr_poly_sin_cos_series(s, c, s, hlen, n, ctx);
            status |= _gr_poly_div_series(f, s, n, c, n, n, ctx);
            status |= gr_neg(t, t, ctx);
            status |= _gr_vec_mul_scalar(f, f, n, t, ctx);
            break;
        /* cot(h) = cos(h)/sin(h) */
        case 2:
            status |= _gr_poly_sin_cos_series(s, c, h, hlen, n, ctx);
            status |= _gr_poly_div_series(f, c, n, s, n, n, ctx);
            break;
        /* coth(h) = i cos(ih)/sin(ih) */
        case 3:
            status |= gr_i(t, ctx);
            status |= _gr_vec_mul_scalar(s, h, hlen, t, ctx);
            status |= _gr_poly_sin_cos_series(s, c, s, hlen, n, ctx);
            status |= _gr_poly_div_series(f, c, n, s, n, n, ctx);
            status |= _gr_vec_mul_scalar(f, f, n, t, ctx);
            break;
        case 4:
            status |= _gr_poly_sin_cos_pi_series(s, c, h, hlen, n, ctx);
            status |= _gr_poly_div_series(f, s, n, c, n, n, ctx);
            break;
        case 6:
            status |= _gr_poly_sin_cos_pi_series(s, c, h, hlen, n, ctx);
            status |= _gr_poly_div_series(f, c, n, s, n, n, ctx);
            break;
    }

    GR_TMP_CLEAR_VEC(s, 2 * n + 1, ctx);

    return status;
}

int
gr_poly_tan_series_sine_cosine(gr_poly_t res, const gr_poly_t h, slong len, int func, gr_ctx_t ctx)
{
    slong hlen = h->length;
    int status = GR_SUCCESS;

    if (hlen == 0 || len == 0)
    {
        if (hlen == 0 && ((func & 3) >= 2))
            return GR_DOMAIN;
        else
            return gr_poly_zero(res, ctx);
    }

    if (hlen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_tan_series_sine_cosine(res->coeffs, h->coeffs, hlen, len, func, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

