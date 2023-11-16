/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

#ifdef __GNUC__
# define log __builtin_log
# define pow __builtin_pow
#else
# include <math.h>
#endif

/* allow changing this from the test code */
FLINT_DLL slong arb_poly_newton_exp_cutoff = 0;

/* with inverse=1 simultaneously computes g = exp(-x) to length n
with inverse=0 uses g as scratch space, computing
g = exp(-x) only to length (n+1)/2 */
static void
_arb_poly_exp_series_newton(arb_ptr f, arb_ptr g,
    arb_srcptr h, slong len, slong prec, int inverse, slong cutoff)
{
    slong alloc;
    arb_ptr T, U, hprime;

    alloc = 3 * len;
    T = _arb_vec_init(alloc);
    U = T + len;
    hprime = U + len;

    _arb_poly_derivative(hprime, h, len, prec);
    arb_zero(hprime + len - 1);

    NEWTON_INIT(cutoff, len)

    /* f := exp(h) + O(x^m), g := exp(-h) + O(x^m2) */
    NEWTON_BASECASE(n)
    _arb_poly_exp_series_basecase(f, h, n, n, prec);
    _arb_poly_inv_series(g, f, (n + 1) / 2, (n + 1) / 2, prec);
    NEWTON_END_BASECASE

    /* extend from length m to length n */
    NEWTON_LOOP(m, n)

    slong m2 = (m + 1) / 2;
    slong l = m - 1; /* shifted for derivative */

    /* g := exp(-h) + O(x^m) */
    _arb_poly_mullow(T, f, m, g, m2, m, prec);
    _arb_poly_mullow(g + m2, g, m2, T + m2, m - m2, m - m2, prec);
    _arb_vec_neg(g + m2, g + m2, m - m2);

    /* U := h' + g (f' - f h') + O(x^(n-1))
        Note: should replace h' by h' mod x^(m-1) */
    _arb_vec_zero(f + m, n - m);
    _arb_poly_mullow(T, f, n, hprime, n, n, prec); /* should be mulmid */
    _arb_poly_derivative(U, f, n, prec); arb_zero(U + n - 1); /* should skip low terms */
    _arb_vec_sub(U + l, U + l, T + l, n - l, prec);
    _arb_poly_mullow(T + l, g, n - m, U + l, n - m, n - m, prec);
    _arb_vec_add(U + l, hprime + l, T + l, n - m, prec);

    /* f := f + f * (h - int U) + O(x^n) = exp(h) + O(x^n) */
    _arb_poly_integral(U, U, n, prec); /* should skip low terms */
    _arb_vec_sub(U + m, h + m, U + m, n - m, prec);
    _arb_poly_mullow(f + m, f, n - m, U + m, n - m, n - m, prec);

    /* g := exp(-h) + O(x^n) */
    /* not needed if we only want exp(x) */
    if (n == len && inverse)
    {
        _arb_poly_mullow(T, f, n, g, m, n, prec);
        _arb_poly_mullow(g + m, g, m, T + m, n - m, n - m, prec);
        _arb_vec_neg(g + m, g + m, n - m);
    }

    NEWTON_END_LOOP

    NEWTON_END

    _arb_vec_clear(T, alloc);
}

void
_arb_poly_exp_series(arb_ptr f, arb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        arb_exp(f, h, prec);
        _arb_vec_zero(f + 1, n - 1);
    }
    else if (n == 2)
    {
        arb_exp(f, h, prec);
        arb_mul(f + 1, f, h + 1, prec);  /* safe since hlen >= 2 */
    }
    else if (_arb_vec_is_zero(h + 1, hlen - 2)) /* h = a + bx^d */
    {
        slong i, j, d = hlen - 1;
        arb_t t;
        arb_init(t);
        arb_set(t, h + d);
        arb_exp(f, h, prec);
        for (i = 1, j = d; j < n; j += d, i++)
        {
            arb_mul(f + j, f + j - d, t, prec);
            arb_div_ui(f + j, f + j, i, prec);
            _arb_vec_zero(f + j - d + 1, hlen - 2);
        }
        _arb_vec_zero(f + j - d + 1, n - (j - d + 1));
        arb_clear(t);
    }
    else
    {
        slong cutoff;

        if (arb_poly_newton_exp_cutoff != 0)
            cutoff = arb_poly_newton_exp_cutoff;
        else if (prec <= 256)
            cutoff = 750;
        else
            cutoff = 1e5 / pow(log(prec), 3);

        if (hlen <= cutoff)
        {
            _arb_poly_exp_series_basecase(f, h, hlen, n, prec);
        }
        else
        {
            arb_ptr g, t;
            arb_t u;
            int fix;

            g = _arb_vec_init((n + 1) / 2);
            fix = (hlen < n || h == f || !arb_is_zero(h));

            if (fix)
            {
                t = _arb_vec_init(n);
                _arb_vec_set(t + 1, h + 1, hlen - 1);
            }
            else
                t = (arb_ptr) h;

            arb_init(u);
            arb_exp(u, h, prec);

            _arb_poly_exp_series_newton(f, g, t, n, prec, 0, cutoff);

            if (!arb_is_one(u))
                _arb_vec_scalar_mul(f, f, n, u, prec);

            _arb_vec_clear(g, (n + 1) / 2);
            if (fix)
                _arb_vec_clear(t, n);
            arb_clear(u);
        }
    }
}

void
arb_poly_exp_series(arb_poly_t f, const arb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        arb_poly_zero(f);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_one(f);
        return;
    }

    if (hlen == 1)
        n = 1;

    arb_poly_fit_length(f, n);
    _arb_poly_exp_series(f->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(f, n);
    _arb_poly_normalise(f);
}
