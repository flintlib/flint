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

/* todo: handle short h */
void
_nmod_poly_atan_series(mp_ptr g, mp_srcptr h, slong n, nmod_t mod)
{
    mp_ptr t, u;

    t = _nmod_vec_init(n);
    u = _nmod_vec_init(n);

    /* atan(h(x)) = integral(h'(x)/(1+h(x)^2)) */
    _nmod_poly_mullow(u, h, n, h, n, n, mod); u[0] = UWORD(1);
    _nmod_poly_derivative(t, h, n, mod); t[n-1] = UWORD(0);
    _nmod_poly_div_series(g, t, n, u, n, n, mod);
    _nmod_poly_integral(g, g, n, mod);

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}

void
nmod_poly_atan_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    mp_ptr h_coeffs;
    slong h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_printf("Exception (nmod_poly_atan_series): Constant term != 0.\n");
        flint_abort();
    }

    if (h_len == 1 || n < 2)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);

    if (h_len < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, h_len);
        flint_mpn_zero(h_coeffs + h_len, n - h_len);
    }
    else
        h_coeffs = h->coeffs;

    _nmod_poly_atan_series(g->coeffs, h_coeffs, n, h->mod);

    if (h_len < n)
        _nmod_vec_clear(h_coeffs);

    g->length = n;
	_nmod_poly_normalise(g);
}
