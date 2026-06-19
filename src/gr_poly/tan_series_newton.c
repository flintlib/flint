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

/*
We use Newton iteration to invert atan; instead of recomputing atan(res)
from scratch each iteration, we do inline Newton iterations to update
the quotient for the arctangent inline. With some sign/scale changes,
we actually support computing the following functions (same as the
basecase algorithm):

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
_gr_poly_tan_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, int func, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, j, m, n;
    gr_ptr Q, u, W, Tprime, pi_val, reciprocals;
    slong a[FLINT_BITS];
    int times_pi, use_reciprocals;

    flen = FLINT_MIN(flen, len);

    times_pi = (func >= 4);

    if (len < cutoff)
        return _gr_poly_tan_series_basecase(res, f, flen, len, func, ctx);

    /* We need m >= 2 throughout the Newton loop so that Tprimelen = m-1 >= 1
       and W2len >= 1.  Raising cutoff to 3 guarantees the basecase precision
       is >= 2, so every m and n in the loop satisfy m >= 2 and n >= 3. */
    cutoff = FLINT_MAX(cutoff, 3);

    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 1) / 2);

    status |= _gr_poly_tan_series_basecase(res, f, flen, n, func, ctx);

    use_reciprocals = (gr_ctx_is_finite_characteristic(ctx) == T_TRUE);

    /* Carry Q = 1/u and u = +/-(1 +/- res^2) alongside res. */
    GR_TMP_INIT_VEC(Q,      len,           ctx);
    GR_TMP_INIT_VEC(u,      len,           ctx);
    GR_TMP_INIT_VEC(Tprime, len - 1,       ctx);
    GR_TMP_INIT_VEC(W,      (len + 1) / 2, ctx);
    if (use_reciprocals)
    {
        GR_TMP_INIT_VEC(reciprocals, len - 1, ctx);
        status |= _gr_vec_reciprocals(reciprocals, len - 1, ctx);
    }

    GR_TMP_INIT(pi_val, ctx);
    if (times_pi)
        status |= gr_pi(pi_val, ctx);


    slong top = i;

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];

        slong Dlen  = FLINT_MIN(2 * m - 1, n);
        slong Wlen  = FLINT_MIN(Dlen + m - 1, n);
        slong W2len = Wlen - m;

        /* Initialise or update u and Q */
        if (i == top - 1)
        {
            /* First Newton step: u = +/-(1 +/- res^2) from scratch. */
            status |= _gr_poly_mullow(u, res, m, res, m, Dlen, ctx);
            if ((func & 3) != 0)
                status |= _gr_vec_neg(u, u, Dlen, ctx);
            status |= gr_add_si(u, u, ((func & 3) == 2) ? -1 : 1, ctx);
            status |= _gr_poly_inv_series(Q, u, Dlen, m, ctx);
        }
        else
        {
            slong m_prev = a[i + 2];
            /* Update high part of u as res changed */
            status |= _gr_poly_mulmid(GR_ENTRY(u, m_prev, sz),
                                      res, m, res, m, m_prev, Dlen, ctx);
            if ((func & 3) != 0)
                status |= _gr_vec_neg(GR_ENTRY(u, m_prev, sz),
                GR_ENTRY(u, m_prev, sz), Dlen - m_prev, ctx);

            /* Extend Q from m_prev to m via a Newton step. */
            slong Wlen_Q  = FLINT_MIN(Dlen + m_prev - 1, m);
            slong W2len_Q = Wlen_Q - m_prev;
            status |= _gr_poly_mulmid(W, u, Dlen, Q, m_prev, m_prev, Wlen_Q, ctx);
            status |= _gr_poly_mullow(GR_ENTRY(Q, m_prev, sz),
                                      Q, m_prev, W, W2len_Q, m - m_prev, ctx);
            status |= _gr_vec_neg(GR_ENTRY(Q, m_prev, sz),
                                  GR_ENTRY(Q, m_prev, sz), m - m_prev, ctx);
        }

        status |= _gr_poly_derivative(Tprime, res, m, ctx);

        /* Extend Q from m to n terms */
        status |= _gr_poly_mulmid(W, u, Dlen, Q, m, m, Wlen, ctx);
        status |= _gr_poly_mullow(GR_ENTRY(Q, m, sz), Q, m, W, W2len, n - m, ctx);
        status |= _gr_vec_neg(GR_ENTRY(Q, m, sz), GR_ENTRY(Q, m, sz), n - m, ctx);

        /* High n-m terms of inv_trig(res) */
        status |= _gr_poly_mulmid(W, Q, n, Tprime, m - 1, m - 1, n - 1, ctx);

        if (use_reciprocals)
            status |= _gr_vec_mul(W, W, GR_ENTRY(reciprocals, m - 1, sz), n - m, ctx);
        else
            for (j = 0; j < n - m; j++)
                status |= gr_div_ui(GR_ENTRY(W, j, sz), GR_ENTRY(W, j, sz), m + j, ctx);

        /* Apply the Newton correction */
        slong fhigh_len = FLINT_MAX(0, FLINT_MIN(flen - m, n - m));
        if (times_pi && fhigh_len > 0)
        {
            status |= _gr_vec_mul_scalar(Tprime, GR_ENTRY(f, m, sz), fhigh_len, pi_val, ctx);
            status |= _gr_poly_sub(Tprime, Tprime, fhigh_len, W, n - m, ctx);
        }
        else
            status |= _gr_poly_sub(Tprime, GR_ENTRY(f, m, sz), fhigh_len, W, n - m, ctx);

        status |= _gr_poly_mullow(GR_ENTRY(res, m, sz), u, Dlen, Tprime, n - m, n - m, ctx);
    }

    GR_TMP_CLEAR(pi_val, ctx);
    GR_TMP_CLEAR_VEC(Q,      len,           ctx);
    GR_TMP_CLEAR_VEC(u,      len,           ctx);
    GR_TMP_CLEAR_VEC(Tprime, len - 1,       ctx);
    GR_TMP_CLEAR_VEC(W,      (len + 1) / 2, ctx);
    if (use_reciprocals)
        GR_TMP_CLEAR_VEC(reciprocals, len - 1, ctx);

    return status;
}

int
gr_poly_tan_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, int func, gr_ctx_t ctx)
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
    status |= _gr_poly_tan_series_newton(res->coeffs, f->coeffs, flen, len, cutoff, func, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

