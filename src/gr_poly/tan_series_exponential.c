/*
    Copyright (C) 2023, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_tan_series_exponential(gr_ptr f, gr_srcptr h, slong hlen, slong n, int func, gr_ctx_t ctx)
{
    gr_ptr e, t, u;
    int status = GR_SUCCESS;
    slong alloc, sz = ctx->sizeof_elem;
    int neg, cot_form, times_pi, is_real_exp;

    hlen = FLINT_MIN(hlen, n);

    cot_form    = ((func & 3) >= 2);
    times_pi    = (func >= 4);
    is_real_exp = (func == 1 || func == 3);

    /*
         tan(h)    = -i * (e^{2ih}  - 1) / (e^{2ih}  + 1)
         tanh(h)   =      (e^{2h}   - 1) / (e^{2h}   + 1)
         cot(h)    = -i * (e^{2ih}  + 1) / (e^{2ih}  - 1)
         coth(h)   =      (e^{2h}   + 1) / (e^{2h}   - 1)
         tan_pi(h) = -i * (e^{2pih} - 1) / (e^{2pih} + 1)
         cot_pi(h) = -i * (e^{2pih} + 1) / (e^{2pih} - 1)
      
      For numerical stability we choose the sign of the exponent so that
      |e^{arg}| < 1, avoiding large exponentials. If csgn or im is not available
      (ambiguous or unsupported), the error is ignored and we proceed with neg = 0.
    */
    {
        gr_ptr sgn;
        GR_TMP_INIT(sgn, ctx);

        if (is_real_exp)
        {
            neg = 1;
            if (gr_csgn(sgn, h, ctx) == GR_SUCCESS)
                if (gr_is_one(sgn, ctx) == T_FALSE)
                    neg = 0;
        }
        else
        {
            neg = 0;
            if (gr_im(sgn, h, ctx) == GR_SUCCESS)
                if (gr_csgn(sgn, sgn, ctx) == GR_SUCCESS)
                    if (gr_is_one(sgn, ctx) == T_FALSE)
                        neg = 1;
        }

        GR_TMP_CLEAR(sgn, ctx);
    }

    alloc = 2 * n + 1;
    GR_TMP_INIT_VEC(e, alloc, ctx);
    t = GR_ENTRY(e, n, sz);
    u = GR_ENTRY(t, n, sz);

    if (is_real_exp)
    {
        status |= gr_set_si(u, neg ? -2 : 2, ctx);
    }
    else if (times_pi)
    {
        status |= gr_pi(u, ctx);
        status |= gr_mul_si(u, u, neg ? -2 : 2, ctx);
        gr_ptr iv;
        GR_TMP_INIT(iv, ctx);
        status |= gr_i(iv, ctx);
        status |= gr_mul(u, u, iv, ctx);
        GR_TMP_CLEAR(iv, ctx);
    }
    else
    {
        status |= gr_i(u, ctx);
        status |= gr_mul_si(u, u, neg ? -2 : 2, ctx);
    }

    status |= _gr_vec_mul_scalar(t, h, hlen, u, ctx);
    status |= _gr_poly_exp_series(e, t, hlen, n, ctx);
    status |= _gr_vec_set(t, e, n, ctx);
    status |= gr_sub_ui(t, t, 1, ctx);
    status |= gr_add_ui(e, e, 1, ctx);

    if (cot_form)
        status |= _gr_poly_div_series(f, e, n, t, n, n, ctx);
    else
        status |= _gr_poly_div_series(f, t, n, e, n, n, ctx);

    if (is_real_exp)
    {
        if (neg)
            status |= _gr_vec_neg(f, f, n, ctx);
    }
    else
    {
        status |= gr_i(u, ctx);
        if (!neg)
            status |= gr_neg(u, u, ctx);
        if (cot_form)
            status |= gr_neg(u, u, ctx);
        status |= _gr_vec_mul_scalar(f, f, n, u, ctx);
    }

    GR_TMP_CLEAR_VEC(e, alloc, ctx);

    return status;
}

int
gr_poly_tan_series_exponential(gr_poly_t res, const gr_poly_t h, slong len, int func, gr_ctx_t ctx)
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
    status |= _gr_poly_tan_series_exponential(res->coeffs, h->coeffs, hlen, len, func, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

