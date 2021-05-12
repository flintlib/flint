/*
    Copyright (C) 2011, 2016, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

#define NMOD_POLY_NEWTON_EXP_CUTOFF 4000

/* c_k x^k -> c_k x^k / (m+k) */
void
_nmod_poly_integral_offset(mp_ptr res, mp_srcptr poly, slong len, slong m, nmod_t mod)
{
    slong k;
    mp_limb_t t;

    t = 1;
    for (k = len - 1; k >= 0; k--)
    {
        res[k] = n_mulmod2_preinv(poly[k], t, mod.n, mod.ninv);
        t = n_mulmod2_preinv(t, m + k, mod.n, mod.ninv);
    }

    if (t >= mod.n)
        t = n_mod2_preinv(t, mod.n, mod.ninv);
    t = n_invmod(t, mod.n);

    for (k = 0; k < len; k++)
    {
        res[k] = n_mulmod2_preinv(res[k], t, mod.n, mod.ninv);
        t = n_mulmod2_preinv(t, m + k, mod.n, mod.ninv);
    }
}

void
_nmod_poly_exp_series_newton(mp_ptr f, mp_ptr g,
    mp_srcptr h, slong hlen, slong n, nmod_t mod)
{
    slong a[FLINT_BITS];
    slong i, m, l, r;
    mp_ptr t, hprime;
    int inverse;

    /* If g is provided, we compute g = exp(-h), and we can use g as
       scratch space. Otherwise, we still need to compute exp(-h) to length
       (n+1)/2 for intermediate use, and we still need n coefficients of
       scratch space. */
    inverse = (g != NULL);
    if (!inverse)
        g = _nmod_vec_init(n);

    hlen = FLINT_MIN(hlen, n);

    t = _nmod_vec_init(n);
    hprime = _nmod_vec_init(hlen - 1);
    _nmod_poly_derivative(hprime, h, hlen, mod);

    for (i = 1; (WORD(1) << i) < n; i++);
    a[i = 0] = n;
    while (n >= NMOD_POLY_NEWTON_EXP_CUTOFF || i == 0)
        a[++i] = (n = (n + 1) / 2);

    /* f := exp(h) + O(x^n),  g := exp(-h) + O(x^n) */
    _nmod_poly_exp_series_basecase(f, h, hlen, n, mod);
    _nmod_poly_inv_series(g, f, n, n, mod);

    for (i--; i >= 0; i--)
    {
        m = n;             /* previous length */
        n = a[i];          /* new length */

        l = FLINT_MIN(hlen, n) - 1;
        r = FLINT_MIN(l + m - 1, n - 1);
        if (l >= m)
            _nmod_poly_mullow(t, hprime, l, f, m, r, mod);
        else
            _nmod_poly_mullow(t, f, m, hprime, l, r, mod);
        _nmod_poly_mullow(g + m, g, n - m, t + m - 1, r + 1 - m, n - m, mod);
        _nmod_poly_integral_offset(g + m, g + m, n - m, m, mod);
        _nmod_poly_mullow(f + m, f, n - m, g + m, n - m, n - m, mod);

        /* g := exp(-h) + O(x^n); not needed if we only want exp(x) */
        if (i != 0 || inverse)
        {
            _nmod_poly_mullow(t, f, n, g, m, n, mod);
            _nmod_poly_mullow(g + m, g, m, t + m, n - m, n - m, mod);
            _nmod_vec_neg(g + m, g + m, n - m, mod);
        }
    }

    _nmod_vec_clear(hprime);
    _nmod_vec_clear(t);
    if (!inverse)
        _nmod_vec_clear(g);
}

void
_nmod_poly_exp_series(mp_ptr f, mp_srcptr h, slong hlen, slong n, nmod_t mod)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen >= 2 && n > 2 && _nmod_vec_is_zero(h + 1, hlen - 2))
        _nmod_poly_exp_series_monomial_ui(f, h[hlen - 1], hlen - 1, n, mod);
    else if (hlen < NMOD_POLY_NEWTON_EXP_CUTOFF)
        _nmod_poly_exp_series_basecase(f, h, hlen, n, mod);
    else
        _nmod_poly_exp_series_newton(f, NULL, h, hlen, n, mod);
}

void 
_nmod_poly_exp_expinv_series(mp_ptr f, mp_ptr g, mp_srcptr h, slong hlen, slong n, nmod_t mod)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen < NMOD_POLY_NEWTON_EXP_CUTOFF)
    {
        _nmod_poly_exp_series_basecase(f, h, hlen, n, mod);
        _nmod_poly_inv_series(g, f, n, n, mod);
    }
    else
    {
        _nmod_poly_exp_series_newton(f, g, h, hlen, n, mod);
    }
}

void
nmod_poly_exp_series(nmod_poly_t f, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_printf("Exception (nmod_poly_exp_series). Constant term != 0.\n");
        flint_abort();
    }

    if (n <= 1 || hlen <= 1)
    {
        if (n == 0)
            nmod_poly_zero(f);
        else
            nmod_poly_one(f);
        return;
    }

    nmod_poly_fit_length(f, n);
    _nmod_poly_exp_series(f->coeffs, h->coeffs, hlen, n, f->mod);
    f->length = n;
    _nmod_poly_normalise(f);
}
