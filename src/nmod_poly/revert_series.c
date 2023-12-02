/*
    Copyright (C) 2010 Sebastian Pancratz
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
#include "gr.h"
#include "gr_poly.h"

void
_nmod_poly_revert_series(mp_ptr Qinv, mp_srcptr Q, slong Qlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_revert_series(Qinv, Q, Qlen, n, ctx));
}

void
nmod_poly_revert_series(nmod_poly_t Qinv,
                                 const nmod_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || Q->coeffs[0] != 0 || Q->coeffs[1] == 0)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_revert_series): "
                "Input must have zero constant and an invertible coefficient of x^1.\n");
    }

    if (Qinv != Q)
    {
        nmod_poly_fit_length(Qinv, n);
        _nmod_poly_revert_series(Qinv->coeffs, Q->coeffs, Q->length, n, Q->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, Q->mod.n, n);
        _nmod_poly_revert_series(t->coeffs, Q->coeffs, Q->length, n, Q->mod);
        nmod_poly_swap(Qinv, t);
        nmod_poly_clear(t);
    }

    _nmod_poly_set_length(Qinv, n);
    _nmod_poly_normalise(Qinv);
}
