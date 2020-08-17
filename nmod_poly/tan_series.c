/*
    Copyright (C) 2011 Fredrik Johansson

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


void
_nmod_poly_tan_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod)
{
    slong m;
    mp_ptr t, u;

    if (n <= 3)
    {
        g[0] = UWORD(0);
        if (n >= 2) g[1] = h[1];
        if (n >= 3) g[2] = h[2];
        return;
    }

    m = (n + 1) / 2;

    _nmod_poly_tan_series(g, h, m, mod);
    _nmod_vec_zero(g + m, n - m);

    t = _nmod_vec_init(n);
    u = _nmod_vec_init(n);

    _nmod_poly_mul(u, g, m, g, m, mod);
    u[0] = UWORD(1);
    if (2*m - 1 < n) u[n-1] = UWORD(0);

    _nmod_poly_atan_series(t, g, n, mod);
    _nmod_vec_sub(t + m, h + m, t + m, n - m, mod);
    _nmod_poly_mullow(g + m, u, n, t + m, n - m, n - m, mod);

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}

void
nmod_poly_tan_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    mp_ptr g_coeffs, h_coeffs;
    nmod_poly_t t1;
    slong h_len;
    
    h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_printf("Exception (nmod_poly_tan_series). Constant term != 0.\n");
        flint_abort();
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

    _nmod_poly_tan_series(g_coeffs, h_coeffs, n, h->mod);

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
