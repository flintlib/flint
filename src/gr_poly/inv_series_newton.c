/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/*
    Newton iteration for 1/Q. With g = 1/Q + O(x^m) and e = 1 - Q g = O(x^m),
    the second-order step is g (1 + e) = 1/Q + O(x^(2m)) and the third-order
    step is g (1 + e + e^2) = 1/Q + O(x^(3m)).

    With FFT multiplication, the third-order step is some 5-10% faster
    for large len. It is only used in finite characteristic: over Z and Q
    the entries of e^2 can be much larger than those of the output (making
    it slower), and in ball arithmetic it can be less numerically stable
    when 1/Q grows faster than Q.
*/
#define INV_NEWTON3_CUTOFF 8192

int
_gr_poly_inv_series_newton(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, m, n, k, Qnlen, Wlen, W2len, alloc;
    gr_ptr W, w, S;
    slong a[FLINT_BITS];
    int order3;

    if (len == 0)
        return GR_SUCCESS;

    if (Qlen == 0)
        return GR_DOMAIN;

    Qlen = FLINT_MIN(Qlen, len);

    if (len < cutoff)
        return _gr_poly_inv_series_basecase(Qinv, Q, Qlen, len, ctx);

    cutoff = FLINT_MAX(cutoff, 2);

    order3 = (len >= INV_NEWTON3_CUTOFF && gr_ctx_is_finite_characteristic(ctx) == T_TRUE);

    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = order3 ? (n + 2) / 3 : (n + 1) / 2);

    status |= _gr_poly_inv_series_basecase(Qinv, Q, Qlen, n, ctx);

    if (status != GR_SUCCESS)
        return status;

    /* Hack: until all rings have a good mulmid, implement both with and without. */
    int have_mulmid = (ctx->methods[GR_METHOD_POLY_MULMID] != (gr_funcptr) _gr_poly_mulmid_generic);

    if (order3)
        alloc = len + 2 * ((len + 2) / 3);
    else
        alloc = have_mulmid ? len / 2 : len;

    GR_TMP_INIT_VEC(W, alloc, ctx);

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];
        k = n - 2 * m;

        Qnlen = FLINT_MIN(Qlen, n);
        Wlen = FLINT_MIN(Qnlen + m - 1, n);
        W2len = Wlen - m;

        FLINT_ASSERT(W2len != 0);
        FLINT_ASSERT(m < Wlen);

        /* w = (Q g)[m, Wlen) = -e / x^m */
        if (have_mulmid)
        {
            status |= _gr_poly_mulmid(W, Q, Qnlen, Qinv, m, m, Wlen, ctx);
            w = W;
        }
        else
        {
            status |= _gr_poly_mullow(W, Q, Qnlen, Qinv, m, Wlen, ctx);
            w = GR_ENTRY(W, m, sz);
        }

        if (k <= 0)
        {
            status |= _gr_poly_mullow(GR_ENTRY(Qinv, m, sz), Qinv, m, w, W2len, n - m, ctx);
            status |= _gr_vec_neg(GR_ENTRY(Qinv, m, sz), GR_ENTRY(Qinv, m, sz), n - m, ctx);
        }
        else
        {
            slong k2, k3;

            /* e + e^2 = -(w - x^m w^2) */
            S = GR_ENTRY(W, have_mulmid ? n - m : n, sz);
            k2 = FLINT_MIN(k, W2len);
            k3 = FLINT_MIN(k, 2 * k2 - 1);
            status |= _gr_poly_mullow(S, w, k2, w, k2, k3, ctx);
            if (k3 < k)
                status |= _gr_vec_zero(GR_ENTRY(S, k3, sz), k - k3, ctx);

            if (W2len >= m)
            {
                /* dense case */
                status |= _gr_vec_zero(GR_ENTRY(w, W2len, sz), n - m - W2len, ctx);
                status |= _gr_vec_sub(GR_ENTRY(w, m, sz), GR_ENTRY(w, m, sz), S, k, ctx);
                status |= _gr_poly_mullow(GR_ENTRY(Qinv, m, sz), Qinv, m, w, n - m, n - m, ctx);
                status |= _gr_vec_neg(GR_ENTRY(Qinv, m, sz), GR_ENTRY(Qinv, m, sz), n - m, ctx);
            }
            else
            {
                /* short Q: both e and e^2 are short */
                slong l1 = FLINT_MIN(n - m, m + W2len - 1);
                slong l2 = FLINT_MIN(k, FLINT_MIN(m, k) + k3 - 1);
                gr_ptr P = GR_ENTRY(S, k, sz);

                status |= _gr_poly_mullow(GR_ENTRY(Qinv, m, sz), Qinv, m, w, W2len, l1, ctx);
                status |= _gr_vec_zero(GR_ENTRY(Qinv, m + l1, sz), n - m - l1, ctx);
                status |= _gr_vec_neg(GR_ENTRY(Qinv, m, sz), GR_ENTRY(Qinv, m, sz), n - m, ctx);
                status |= _gr_poly_mullow(P, Qinv, FLINT_MIN(m, k), S, k3, l2, ctx);
                status |= _gr_vec_add(GR_ENTRY(Qinv, 2 * m, sz), GR_ENTRY(Qinv, 2 * m, sz), P, l2, ctx);
            }
        }
    }

    GR_TMP_CLEAR_VEC(W, alloc, ctx);

    return status;
}

int
gr_poly_inv_series_newton(gr_poly_t Qinv, const gr_poly_t Q, slong len, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong Qlen;

    if (len == 0)
        return gr_poly_zero(Qinv, ctx);

    Qlen = Q->length;

    if (Qlen == 0)
        return GR_DOMAIN;

    if (Qlen == 1)
        len = 1;

    if (Qinv == Q)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_inv_series_newton(t, Q, len, cutoff, ctx);
        gr_poly_swap(Qinv, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Qinv, len, ctx);
    status |= _gr_poly_inv_series_newton(Qinv->coeffs, Q->coeffs, Q->length, len, cutoff, ctx);
    _gr_poly_set_length_normalise(Qinv, len, ctx);
    return status;
}
