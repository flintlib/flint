/*
    Copyright (C) 2011, 2016, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_exp_series(mp_ptr f, mp_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_exp_series(f, h, hlen, n, ctx));
}

/* todo: gr version */
void
_nmod_poly_exp_expinv_series(mp_ptr f, mp_ptr g, mp_srcptr h, slong hlen, slong n, nmod_t mod)
{
    _nmod_poly_exp_series(f, h, hlen, n, mod);
    _nmod_poly_inv_series(g, f, n, n, mod);
}

void
nmod_poly_exp_series(nmod_poly_t f, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_exp_series). Constant term != 0.\n");
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
