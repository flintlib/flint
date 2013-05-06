/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

#define NMOD_POLY_NEWTON_EXP_CUTOFF 200

/* with inverse=1 simultaneously computes g = exp(-x) to length n
   with inverse=0 uses g as scratch space, computing
                  g = exp(-x) only to length (n+1)/2 */
static void
_nmod_poly_exp_series_newton(mp_ptr f, mp_ptr g,
    mp_srcptr h, long n, nmod_t mod, int inverse)
{
    long a[FLINT_BITS];
    long i, m, m2, l;
    mp_ptr T, U, hprime;

    T = _nmod_vec_init(n);
    U = _nmod_vec_init(n);
    hprime = _nmod_vec_init(n);

    _nmod_poly_derivative(hprime, h, n, mod);
    hprime[n-1] = 0UL;

    for (i = 1; (1L << i) < n; i++);
    a[i = 0] = n;
    while (n >= NMOD_POLY_NEWTON_EXP_CUTOFF)
        a[++i] = (n = (n + 1) / 2);

    /* f := exp(h) + O(x^m),  g := exp(-h) + O(x^m2) */
    _nmod_poly_exp_series_basecase(f, h, n, n, mod);
    _nmod_poly_inv_series_basecase(g, f, (n + 1) / 2, mod);

    for (i--; i >= 0; i--)
    {
        m = n;             /* previous length */
        n = a[i];          /* new length */
        m2 = (m + 1) / 2;
        l = m - 1;         /* shifted for derivative */

        /* g := exp(-h) + O(x^m) */
        _nmod_poly_mullow(T, f, m, g, m2, m, mod);
        _nmod_poly_mullow(g + m2, g, m2, T + m2, m - m2, m - m2, mod);
        _nmod_vec_neg(g + m2, g + m2, m - m2, mod);

        /* U := h' + g (f' - f h') + O(x^(n-1))
           Note: should replace h' by h' mod x^(m-1) */
        _nmod_vec_zero(f + m, n - m);
        _nmod_poly_mullow(T, f, n, hprime, n, n, mod);     /* should be mulmid */
        _nmod_poly_derivative(U, f, n, mod); U[n - 1] = 0; /* should skip low terms */
        _nmod_vec_sub(U + l, U + l, T + l, n - l, mod);
        _nmod_poly_mullow(T + l, g, n - m, U + l, n - m, n - m, mod);
        _nmod_vec_add(U + l, hprime + l, T + l, n - m, mod);

        /* f := f + f * (h - int U) + O(x^n) = exp(h) + O(x^n) */
        _nmod_poly_integral(U, U, n, mod);  /* should skip low terms */
        _nmod_vec_sub(U + m, h + m, U + m, n - m, mod);
        _nmod_poly_mullow(f + m, f, n - m, U + m, n - m, n - m, mod);

        /* g := exp(-h) + O(x^n) */
        /* not needed if we only want exp(x) */
        if (i == 0 && inverse)
        {
            _nmod_poly_mullow(T, f, n, g, m, n, mod);
            _nmod_poly_mullow(g + m, g, m, T + m, n - m, n - m, mod);
            _nmod_vec_neg(g + m, g + m, n - m, mod);
        }
    }

    _nmod_vec_clear(hprime);
    _nmod_vec_clear(T);
    _nmod_vec_clear(U);
}

void 
_nmod_poly_exp_expinv_series(mp_ptr f, mp_ptr g, mp_srcptr h, long n, nmod_t mod)
{
    if (n < NMOD_POLY_NEWTON_EXP_CUTOFF)
    {
        _nmod_poly_exp_series_basecase(f, h, n, n, mod);
        _nmod_poly_inv_series(g, f, n, mod);
    }
    else
    {
        _nmod_poly_exp_series_newton(f, g, h, n, mod, 1);
    }
}

void 
_nmod_poly_exp_series(mp_ptr f, mp_srcptr h, long n, nmod_t mod)
{
    if (n < NMOD_POLY_NEWTON_EXP_CUTOFF)
    {
        _nmod_poly_exp_series_basecase(f, h, n, n, mod);
    }
    else
    {
        mp_ptr g = _nmod_vec_init((n + 1) / 2);
        _nmod_poly_exp_series_newton(f, g, h, n, mod, 0);
        _nmod_vec_clear(g);
    }
}

void
nmod_poly_exp_series(nmod_poly_t f, const nmod_poly_t h, long n)
{
    mp_ptr f_coeffs, h_coeffs;
    nmod_poly_t t1;
    long hlen, k;

    nmod_poly_fit_length(f, n);
    hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != 0UL)
    {
        printf("Exception (nmod_poly_exp_series). Constant term != 0.\n");
        abort();
    }

    if (n <= 1 || hlen == 0)
    {
        if (n == 0)
        {
            nmod_poly_zero(f);
        }
        else
        {
            f->coeffs[0] = 1UL;
            f->length = 1;
        }
        return;
    }

    /* Handle monomials */
    for (k = 0; h->coeffs[k] == 0UL && k < n - 1; k++);
    if (k == hlen - 1 || k == n - 1)
    {
        hlen = FLINT_MIN(hlen, n);
        _nmod_poly_exp_series_monomial_ui(f->coeffs,
            h->coeffs[hlen-1], hlen - 1, n, f->mod);
        f->length = n;
        _nmod_poly_normalise(f);
        return;
    }

    if (n < NMOD_POLY_NEWTON_EXP_CUTOFF)
    {
        _nmod_poly_exp_series_basecase(f->coeffs, h->coeffs, hlen, n, f->mod);
        f->length = n;
        _nmod_poly_normalise(f);
        return;
    }

    if (hlen < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, hlen);
        flint_mpn_zero(h_coeffs + hlen, n - hlen);
    }
    else
        h_coeffs = h->coeffs;

    if (h == f && hlen >= n)
    {
        nmod_poly_init2(t1, h->mod.n, n);
        f_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(f, n);
        f_coeffs = f->coeffs;
    }

    _nmod_poly_exp_series(f_coeffs, h_coeffs, n, f->mod);

    if (h == f && hlen >= n)
    {
        nmod_poly_swap(f, t1);
        nmod_poly_clear(t1);
    }
    
    f->length = n;

    if (hlen < n)
        _nmod_vec_clear(h_coeffs);

    _nmod_poly_normalise(f);
}
