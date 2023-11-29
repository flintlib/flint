/*
    Copyright (C) 2011, 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_log_series(mp_ptr g, mp_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_log_series(g, h, hlen, n, ctx));
}

void
nmod_poly_log_series(nmod_poly_t res, const nmod_poly_t f, slong n)
{
    slong flen = f->length;

    if (flen < 1 || f->coeffs[0] != UWORD(1))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_log_series). Constant term != 1.\n");
    }

    if (flen == 1 || n < 2)
    {
        nmod_poly_zero(res);
        return;
    }

    nmod_poly_fit_length(res, n);
    _nmod_poly_log_series(res->coeffs, f->coeffs, f->length, n, res->mod);
    res->length = n;
	_nmod_poly_normalise(res);
}

