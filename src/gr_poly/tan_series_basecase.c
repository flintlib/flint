/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"

/*
The Taylor coefficients of tan(h) satisfy the triangular recurrence

    t_k = h_k (1 + t_0^2) + (1/k) \sum_{j=1}^{k-1} j h_j C_{k-j}
    C_m = \sum_{i=0}^{m} t_i * t_{m-i}

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
_gr_poly_tan_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, int func, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong k, sz = ctx->sizeof_elem;
    gr_ptr tmp, v, C, hprime, h_pi, reciprocals;
    gr_srcptr h_or_h_pi;
    slong alloc;
    int use_reciprocals, negate, times_pi;

    hlen = FLINT_MIN(hlen, n);

    switch (func)
    {
        case 0: status |= gr_tan(f, h, ctx); break;
        case 1: status |= gr_tanh(f, h, ctx); break;
        case 2: status |= gr_cot(f, h, ctx); break;
        case 3: status |= gr_coth(f, h, ctx); break;
        case 4: status |= gr_tan_pi(f, h, ctx); break;
        case 6: status |= gr_cot_pi(f, h, ctx); break;
    }

    if (status != GR_SUCCESS)
        return status;

    if (hlen == 1)
        return _gr_vec_zero(GR_ENTRY(f, 1, sz), n - 1, ctx);

    negate = ((func & 3) >= 1);
    times_pi = (func >= 4);

    if (n == 2)
    {
        GR_TMP_INIT(v, ctx);
        status |= gr_mul(v, f, f, ctx);
        status |= gr_add_si(v, v, (func & 1) ? -1 : 1, ctx);
        if (negate)
            status |= gr_neg(v, v, ctx);
        status |= gr_mul(GR_ENTRY(f, 1, sz), GR_ENTRY(h, 1, sz), v, ctx);
        if (times_pi)
        {
            status |= gr_pi(v, ctx);
            status |= gr_mul(GR_ENTRY(f, 1, sz), GR_ENTRY(f, 1, sz), v, ctx);
        }
        GR_TMP_CLEAR(v, ctx);
        return status;
    }

    use_reciprocals = (gr_ctx_is_finite_characteristic(ctx) == T_TRUE);

    alloc = (hlen - 1) + (n) + (1) + (times_pi ? hlen : 0) + (use_reciprocals ? (n - 1) : 0);
    GR_TMP_INIT_VEC(tmp, alloc, ctx);

    hprime = tmp;
    C = GR_ENTRY(hprime, hlen - 1, sz);
    v = GR_ENTRY(C, n, sz);
    h_pi = GR_ENTRY(v, 1, sz);
    reciprocals = GR_ENTRY(h_pi, (times_pi ? hlen : 0), sz);

    if (times_pi)
    {
        status |= gr_pi(h_pi, ctx);
        status |= _gr_vec_mul_scalar(GR_ENTRY(h_pi, 1, sz), GR_ENTRY(h, 1, sz), hlen - 1, h_pi, ctx);
        h_or_h_pi = h_pi;
    }
    else
    {
        h_or_h_pi = h;
    }

    status |= gr_mul(C, f, f, ctx);
    status |= gr_add_si(v, C, (func & 1) ? -1 : 1, ctx);
    if (negate)
        status |= gr_neg(v, v, ctx);

    status |= _gr_poly_derivative(hprime, h_or_h_pi, hlen, ctx);

    if (use_reciprocals)
        status |= _gr_vec_reciprocals(reciprocals, n - 1, ctx);

    /* Unroll k = 1. */
    status |= gr_mul(GR_ENTRY(f, 1, sz), GR_ENTRY(h_or_h_pi, 1, sz), v, ctx);
    if (n > 2)
    {
        status |= gr_mul(GR_ENTRY(C, 1, sz), GR_ENTRY(f, 0, sz), GR_ENTRY(f, 1, sz), ctx);
        status |= gr_mul_two(GR_ENTRY(C, 1, sz), GR_ENTRY(C, 1, sz), ctx);
    }

    for (k = 2; k < n && status == GR_SUCCESS; k++)
    {
        slong m = FLINT_MIN(k - 1, hlen - 1);

        gr_ptr Ck = GR_ENTRY(C, k, sz);
        gr_ptr fk = GR_ENTRY(f, k, sz);
        /* Ck is available as scratch space before we update it below */
        gr_ptr t = Ck;

        if (k < hlen)
            status |= gr_mul(t, GR_ENTRY(h_or_h_pi, k, sz), v, ctx);

        status |= _gr_vec_dot_rev(fk, NULL, negate, hprime, GR_ENTRY(C, k - m, sz), m, ctx);

        if (use_reciprocals)
            status |= gr_mul(fk, fk, GR_ENTRY(reciprocals, k - 1, sz), ctx);
        else
            status |= gr_div_ui(fk, fk, k, ctx);

        if (k < hlen)
            status |= gr_add(fk, fk, t, ctx);

        /* update Ck for next iteration */
        if (k < n - 1)
        {
            slong half = (k + 1) / 2;
            status |= _gr_vec_dot_rev(Ck, NULL, 0,
                                      f, GR_ENTRY(f, k - half + 1, sz), half, ctx);
            status |= gr_mul_two(Ck, Ck, ctx);
            if (k % 2 == 0)
                status |= gr_addmul(GR_ENTRY(C, k, sz),
                                    GR_ENTRY(f, half, sz), GR_ENTRY(f, half, sz), ctx);
        }
    }

    GR_TMP_CLEAR_VEC(tmp, alloc, ctx);

    return status;
}

int
gr_poly_tan_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, int func, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
    {
        if (f->length == 0 && ((func & 3) >= 2))
            return GR_DOMAIN;
        else
            return gr_poly_zero(res, ctx);
    }

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_tan_series_basecase(res->coeffs, f->coeffs, flen, len, func, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

