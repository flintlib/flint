/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_sinh_series(mp_ptr f, mp_srcptr h, slong n, nmod_t mod)
{
    mp_ptr g = _nmod_vec_init(n);
    _nmod_poly_exp_expinv_series(f, g, h, n, n, mod);
    _nmod_vec_sub(f, f, g, n, mod);
    _nmod_vec_scalar_mul_nmod(f, f, n, n_invmod(UWORD(2), mod.n), mod);
    _nmod_vec_clear(g);
}

void
nmod_poly_sinh_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    mp_ptr g_coeffs, h_coeffs;
    nmod_poly_t t1;
    slong h_len;

    h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_sinh_series). Constant term != 0.\n");
    }

    if (h_len == 1 || n < 2)
    {
        nmod_poly_zero(g);
        return;
    }

    if (h_len < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, h_len);
        flint_mpn_zero(h_coeffs + h_len, n - h_len);
    }
    else
        h_coeffs = h->coeffs;

    if (h == g && h_len >= n)
    {
        nmod_poly_init2(t1, h->mod.n, n);
        g_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(g, n);
        g_coeffs = g->coeffs;
    }

    _nmod_poly_sinh_series(g_coeffs, h_coeffs, n, h->mod);

    if (h == g && h_len >= n)
    {
        nmod_poly_swap(g, t1);
        nmod_poly_clear(t1);
    }

    g->length = n;

    if (h_len < n)
        _nmod_vec_clear(h_coeffs);

    _nmod_poly_normalise(g);
}
