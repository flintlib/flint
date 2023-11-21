/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

/* (a + bx^c)^g where a = f[0] and b = f[flen-1] */
void
_acb_poly_binomial_pow_acb_series(acb_ptr h, acb_srcptr f, slong flen, const acb_t g, slong len, slong prec)
{
    slong i, j, d;
    acb_t t;

    acb_init(t);

    d = flen - 1;
    acb_pow(h, f, g, prec);
    acb_div(t, f + d, f, prec);

    for (i = 1, j = d; j < len; i++, j += d)
    {
        acb_sub_ui(h + j, g, i - 1, prec);
        acb_mul(h + j, h + j, h + j - d, prec);
        acb_mul(h + j, h + j, t, prec);
        acb_div_ui(h + j, h + j, i, prec);
    }

    if (d > 1)
    {
        for (i = 1; i < len; i++)
            if (i % d != 0)
                acb_zero(h + i);
    }

    acb_clear(t);
    return;
}

void
_acb_poly_pow_acb_series(acb_ptr h,
    acb_srcptr f, slong flen, const acb_t g, slong len, slong prec)
{
    int f_binomial, g_exact, g_int;

    while (flen > 0 && acb_is_zero(f + flen - 1))
        flen--;

    if (flen <= 1)
    {
        acb_pow(h, f, g, prec);
        _acb_vec_zero(h + 1, len - 1);
        return;
    }

    g_exact = acb_is_exact(g);
    g_int = acb_is_real(g) && arb_is_int(acb_realref(g));
    f_binomial = _acb_vec_is_zero(f + 1, flen - 2);

    /* g = small integer */
    if (g_exact && g_int &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(g)), FLINT_BITS - 1) < 0)
    {
        slong e, hlen;

        e = arf_get_si(arb_midref(acb_realref(g)), ARF_RND_DOWN);
        hlen = poly_pow_length(flen, FLINT_ABS(e), len);

        if (e >= 0)
        {
            _acb_poly_pow_ui_trunc_binexp(h, f, flen, e, hlen, prec);
            _acb_vec_zero(h + hlen, len - hlen);
            return;
        }
        else if (!f_binomial)
        {
            acb_ptr t;
            t = _acb_vec_init(hlen);
            _acb_poly_pow_ui_trunc_binexp(t, f, flen, -e, hlen, prec);
            _acb_poly_inv_series(h, t, hlen, len, prec);
            _acb_vec_clear(t, hlen);
            return;
        }
    }

    /* (a + bx^c)^g */
    if (f_binomial)
    {
        _acb_poly_binomial_pow_acb_series(h, f, flen, g, len, prec);
        return;
    }

    /* g = +/- 1/2 */
    if (g_exact && acb_is_real(g) && arf_cmpabs_2exp_si(arb_midref(acb_realref(g)), -1) == 0)
    {
        if (arf_sgn(arb_midref(acb_realref(g))) > 0)
            _acb_poly_sqrt_series(h, f, flen, len, prec);
        else
            _acb_poly_rsqrt_series(h, f, flen, len, prec);
        return;
    }

    /* f^g = exp(g*log(f)) */
    _acb_poly_log_series(h, f, flen, len, prec);
    _acb_vec_scalar_mul(h, h, len, g, prec);
    _acb_poly_exp_series(h, h, len, len, prec);

}

void
acb_poly_pow_acb_series(acb_poly_t h,
    const acb_poly_t f, const acb_t g, slong len, slong prec)
{
    slong flen;

    flen = f->length;
    flen = FLINT_MIN(flen, len);

    if (len == 0)
    {
        acb_poly_zero(h);
        return;
    }

    if (acb_is_zero(g))
    {
        acb_poly_one(h);
        return;
    }

    if (flen == 0)
    {
        acb_poly_zero(h);
        return;
    }

    if (f == h)
    {
        acb_poly_t t;
        acb_poly_init2(t, len);
        _acb_poly_pow_acb_series(t->coeffs, f->coeffs, flen, g, len, prec);
        _acb_poly_set_length(t, len);
        _acb_poly_normalise(t);
        acb_poly_swap(t, h);
        acb_poly_clear(t);
    }
    else
    {
        acb_poly_fit_length(h, len);
        _acb_poly_pow_acb_series(h->coeffs, f->coeffs, flen, g, len, prec);
        _acb_poly_set_length(h, len);
        _acb_poly_normalise(h);
    }
}

