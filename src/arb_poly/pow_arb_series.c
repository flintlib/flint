/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

/* (a + bx^c)^g where a = f[0] and b = f[flen-1] */
void
_arb_poly_binomial_pow_arb_series(arb_ptr h, arb_srcptr f, slong flen, const arb_t g, slong len, slong prec)
{
    slong i, j, d;
    arb_t t;

    arb_init(t);

    d = flen - 1;
    arb_pow(h, f, g, prec);
    arb_div(t, f + d, f, prec);

    for (i = 1, j = d; j < len; i++, j += d)
    {
        arb_sub_ui(h + j, g, i - 1, prec);
        arb_mul(h + j, h + j, h + j - d, prec);
        arb_mul(h + j, h + j, t, prec);
        arb_div_ui(h + j, h + j, i, prec);
    }

    if (d > 1)
    {
        for (i = 1; i < len; i++)
            if (i % d != 0)
                arb_zero(h + i);
    }

    arb_clear(t);
    return;
}

void
_arb_poly_pow_arb_series(arb_ptr h,
    arb_srcptr f, slong flen, const arb_t g, slong len, slong prec)
{
    int f_binomial, g_exact, g_int;

    while (flen > 0 && arb_is_zero(f + flen - 1))
        flen--;

    if (flen <= 1)
    {
        arb_pow(h, f, g, prec);
        _arb_vec_zero(h + 1, len - 1);
        return;
    }

    g_exact = arb_is_exact(g);
    g_int = arb_is_int(g);
    f_binomial = _arb_vec_is_zero(f + 1, flen - 2);

    /* g = small integer */
    if (g_exact && g_int &&
            arf_cmpabs_2exp_si(arb_midref(g), FLINT_BITS - 1) < 0)
    {
        slong e, hlen;

        e = arf_get_si(arb_midref(g), ARF_RND_DOWN);
        hlen = poly_pow_length(flen, FLINT_ABS(e), len);

        if (e >= 0)
        {
            _arb_poly_pow_ui_trunc_binexp(h, f, flen, e, hlen, prec);
            _arb_vec_zero(h + hlen, len - hlen);
            return;
        }
        else if (!f_binomial)
        {
            arb_ptr t;
            t = _arb_vec_init(hlen);
            _arb_poly_pow_ui_trunc_binexp(t, f, flen, -e, hlen, prec);
            _arb_poly_inv_series(h, t, hlen, len, prec);
            _arb_vec_clear(t, hlen);
            return;
        }
    }

    /* (a + bx^c)^g */


    if (f_binomial)
    {
        _arb_poly_binomial_pow_arb_series(h, f, flen, g, len, prec);
        return;
    }

    /* g = +/- 1/2 */
    if (g_exact && arf_cmpabs_2exp_si(arb_midref(g), -1) == 0)
    {
        if (arf_sgn(arb_midref(g)) > 0)
            _arb_poly_sqrt_series(h, f, flen, len, prec);
        else
            _arb_poly_rsqrt_series(h, f, flen, len, prec);
        return;
    }

    /* f^g = exp(g*log(f)) */
    _arb_poly_log_series(h, f, flen, len, prec);
    _arb_vec_scalar_mul(h, h, len, g, prec);
    _arb_poly_exp_series(h, h, len, len, prec);

}

void
arb_poly_pow_arb_series(arb_poly_t h,
    const arb_poly_t f, const arb_t g, slong len, slong prec)
{
    slong flen;

    flen = f->length;
    flen = FLINT_MIN(flen, len);

    if (len == 0)
    {
        arb_poly_zero(h);
        return;
    }

    if (arb_is_zero(g))
    {
        arb_poly_one(h);
        return;
    }

    if (flen == 0)
    {
        arb_poly_zero(h);
        return;
    }

    if (f == h)
    {
        arb_poly_t t;
        arb_poly_init2(t, len);
        _arb_poly_pow_arb_series(t->coeffs, f->coeffs, flen, g, len, prec);
        _arb_poly_set_length(t, len);
        _arb_poly_normalise(t);
        arb_poly_swap(t, h);
        arb_poly_clear(t);
    }
    else
    {
        arb_poly_fit_length(h, len);
        _arb_poly_pow_arb_series(h->coeffs, f->coeffs, flen, g, len, prec);
        _arb_poly_set_length(h, len);
        _arb_poly_normalise(h);
    }
}

