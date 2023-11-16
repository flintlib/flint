/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"
#include "ca_poly.h"

ca_field_ptr
_ca_vec_same_field2(ca_srcptr A, slong Alen, ca_srcptr B, slong Blen, ca_ctx_t ctx);

void
_ca_poly_exp_series_basecase(ca_ptr f,
        ca_srcptr h, slong hlen, slong len, ca_ctx_t ctx)
{
    slong k;
    ca_ptr a;
    ca_t s, e;

    hlen = FLINT_MIN(hlen, len);

    ca_init(e, ctx);
    ca_exp(e, h, ctx);

    if (_ca_vec_is_fmpq_vec(h + 1, hlen - 1, ctx))
    {
        fmpz *p, *r;
        fmpz_t pden, rden;

        p = _fmpz_vec_init(hlen);
        r = _fmpz_vec_init(len);
        fmpz_init(pden);
        fmpz_init(rden);

        _ca_vec_fmpq_vec_get_fmpz_vec_den(p + 1, pden, h + 1, hlen - 1, ctx);
        _fmpq_poly_exp_series(r, rden, p, pden, hlen, len);
        _ca_vec_set_fmpz_vec_div_fmpz(f, r, rden, len, ctx);

        fmpz_clear(pden);
        fmpz_clear(rden);
        _fmpz_vec_clear(p, hlen);
        _fmpz_vec_clear(r, len);
    }
    else
    {
        ca_init(s, ctx);
        a = _ca_vec_init(hlen, ctx);

        for (k = 1; k < hlen; k++)
            ca_mul_ui(a + k, h + k, k, ctx);

        ca_one(f, ctx);

        for (k = 1; k < len; k++)
        {
            ca_dot(s, NULL, 0, a + 1, 1, f + k - 1, -1, FLINT_MIN(k, hlen - 1), ctx);
            ca_div_ui(f + k, s, k, ctx);
        }

        _ca_vec_clear(a, hlen, ctx);
        ca_clear(s, ctx);
    }

    ca_swap(f, e, ctx);
    _ca_vec_scalar_mul_ca(f + 1, f + 1, len - 1, f, ctx);

    ca_clear(e, ctx);
}

/* c_k x^k -> c_k x^k / (m+k) */
static void
_ca_poly_integral_offset(ca_ptr res, ca_srcptr poly, slong len, slong m, ca_ctx_t ctx)
{
    slong k;

    for (k = 0; k < len; k++)
        ca_div_ui(res + k, poly + k, m + k, ctx);
}

void
_ca_poly_exp_series_newton(ca_ptr f, ca_ptr g,
    ca_srcptr h, slong hlen, slong n, ca_ctx_t ctx)
{
    slong a[FLINT_BITS];
    slong i, m, l, r, alloc;
    ca_ptr t, hprime;
    int inverse;

    if (!(CA_IS_QQ(h, ctx) && fmpq_is_zero(CA_FMPQ(h))))
    {
        hlen = FLINT_MIN(hlen, n);
        t = _ca_vec_init(hlen + 1, ctx);
        ca_exp(t + hlen, h, ctx);
        _ca_vec_set(t + 1, h + 1, hlen - 1, ctx);
        _ca_poly_exp_series_newton(f, g, t, hlen, n, ctx);
        _ca_vec_scalar_mul_ca(f, f, n, t + hlen, ctx);
        if (g != NULL)
            _ca_vec_scalar_div_ca(g, g, n, t + hlen, ctx);
        _ca_vec_clear(t, hlen + 1, ctx);
        return;
    }

    /* If g is provided, we compute g = exp(-h), and we can use g as
       scratch space. Otherwise, we still need to compute exp(-h) to length
       (n+1)/2 for intermediate use, and we still need n coefficients of
       scratch space. */
    alloc = n;
    inverse = (g != NULL);
    if (!inverse)
        g = _ca_vec_init(n, ctx);

    hlen = FLINT_MIN(hlen, n);

    t = _ca_vec_init(n, ctx);
    hprime = _ca_vec_init(hlen - 1, ctx);
    _ca_poly_derivative(hprime, h, hlen, ctx);

    for (i = 1; (WORD(1) << i) < n; i++);
    a[i = 0] = n;
    while (n >= 15 || i == 0)
        a[++i] = (n = (n + 1) / 2);

    /* f := exp(h) + O(x^n),  g := exp(-h) + O(x^n) */
    _ca_poly_exp_series_basecase(f, h, FLINT_MIN(hlen, n), n, ctx);
    _ca_poly_inv_series(g, f, n, n, ctx);

    for (i--; i >= 0; i--)
    {
        m = n;             /* previous length */
        n = a[i];          /* new length */

        l = FLINT_MIN(hlen, n) - 1;
        r = FLINT_MIN(l + m - 1, n - 1);
        if (l >= m)
            _ca_poly_mullow(t, hprime, l, f, m, r, ctx);
        else
            _ca_poly_mullow(t, f, m, hprime, l, r, ctx);
        _ca_poly_mullow(g + m, g, n - m, t + m - 1, r + 1 - m, n - m, ctx);
        _ca_poly_integral_offset(g + m, g + m, n - m, m, ctx);
        _ca_poly_mullow(f + m, f, n - m, g + m, n - m, n - m, ctx);

        /* g := exp(-h) + O(x^n); not needed if we only want exp(x) */
        if (i != 0 || inverse)
        {
            _ca_poly_mullow(t, f, n, g, m, n, ctx);
            _ca_poly_mullow(g + m, g, m, t + m, n - m, n - m, ctx);
            _ca_vec_neg(g + m, g + m, n - m, ctx);
        }
    }

    _ca_vec_clear(hprime, hlen - 1, ctx);
    _ca_vec_clear(t, alloc, ctx);
    if (!inverse)
        _ca_vec_clear(g, alloc, ctx);
}

void
_ca_poly_exp_series(ca_ptr f, ca_srcptr h, slong hlen, slong len, ca_ctx_t ctx)
{
    hlen = FLINT_MIN(hlen, len);

    if (CA_IS_SPECIAL(h))
    {
        if (ca_is_unknown(h, ctx))
            _ca_vec_unknown(f, len, ctx);
        else
            _ca_vec_undefined(f, len, ctx);
        return;
    }

    if (hlen == 1)
    {
        ca_exp(f, h, ctx);
        _ca_vec_zero(f + 1, len - 1, ctx);
    }
    else if (len == 2)
    {
        ca_exp(f, h, ctx);
        ca_mul(f + 1, f, h + 1, ctx);  /* safe since hlen >= 2 */
    }
    else if (_ca_vec_check_is_zero(h + 1, hlen - 2, ctx) == T_TRUE) /* h = a + bx^d */
    {
        slong i, j, d = hlen - 1;
        ca_t t;
        ca_init(t, ctx);
        ca_set(t, h + d, ctx);
        ca_exp(f, h, ctx);
        for (i = 1, j = d; j < len; j += d, i++)
        {
            ca_mul(f + j, f + j - d, t, ctx);
            ca_div_ui(f + j, f + j, i, ctx);
            _ca_vec_zero(f + j - d + 1, hlen - 2, ctx);
        }
        _ca_vec_zero(f + j - d + 1, len - (j - d + 1), ctx);
        ca_clear(t, ctx);
    }
    else
    {
        if (hlen >= 8)
        {
            ca_field_ptr K;

            K = _ca_vec_same_field2(h + 1, hlen - 1, NULL, 0, ctx);

            /* Newton iteration where we have fast multiplication */
            if (K != NULL && CA_FIELD_IS_NF(K))
            {
                if (len >= qqbar_degree(CA_FIELD_NF_QQBAR(K)))
                {
                    _ca_poly_exp_series_newton(f, NULL, h, hlen, len, ctx);
                    return;
                }
            }
        }

        _ca_poly_exp_series_basecase(f, h, hlen, len, ctx);
    }
}

void
ca_poly_exp_series(ca_poly_t f, const ca_poly_t h, slong len, ca_ctx_t ctx)
{
    slong hlen = h->length;

    if (len == 0)
    {
        ca_poly_zero(f, ctx);
        return;
    }

    if (hlen == 0)
    {
        ca_poly_one(f, ctx);
        return;
    }

    if (hlen == 1 && ca_check_is_number(h->coeffs, ctx) == T_TRUE)
        len = 1;

    ca_poly_fit_length(f, len, ctx);
    _ca_poly_exp_series(f->coeffs, h->coeffs, hlen, len, ctx);
    _ca_poly_set_length(f, len, ctx);
    _ca_poly_normalise(f, ctx);
}
