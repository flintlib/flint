/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_atan_series(mp_ptr g, mp_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_atan_series(g, h, hlen, n, ctx));
}

void
nmod_poly_atan_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_atan_series). Constant term != 0.\n");
    }

    if (hlen <= 1 || n <= 1)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);
    _nmod_poly_atan_series(g->coeffs, h->coeffs, hlen, n, h->mod);
    _nmod_poly_set_length(g, n);
	_nmod_poly_normalise(g);
}
