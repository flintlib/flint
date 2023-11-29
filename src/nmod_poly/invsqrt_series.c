/*
    Copyright (C) 2011 William Hart
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
_nmod_poly_invsqrt_series(mp_ptr g, mp_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_rsqrt_series(g, h, hlen, n, ctx));
}

void
nmod_poly_invsqrt_series(nmod_poly_t res, const nmod_poly_t h, slong len)
{
    slong hlen;
    hlen = h->length;

    if (h->length == 0 || h->coeffs[0] == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_invsqrt_series). Division by zero.\n");
    }

    if (len == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (hlen == 1)
        len = 1;

    if (res == h)
    {
        nmod_poly_t t;
        nmod_poly_init_preinv(t, h->mod.n, h->mod.ninv);
        nmod_poly_invsqrt_series(t, h, len);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }
    else
    {
        nmod_poly_fit_length(res, len);
        _nmod_poly_invsqrt_series(res->coeffs, h->coeffs, h->length, len, h->mod);
        _nmod_poly_set_length(res, len);
        _nmod_poly_normalise(res);
    }
}
