/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_hypgeom.h"
#include "gr_poly.h"

void
_arb_hypgeom_bessel_j_series(arb_ptr res, const arb_t nu,
    arb_srcptr z, slong zlen, slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    if (_gr_poly_bessel_j_series(res, nu, z, zlen, len, ctx) != GR_SUCCESS)
        _arb_vec_indeterminate(res, len);
}

void
arb_hypgeom_bessel_j_series(arb_poly_t res, const arb_t nu,
    const arb_poly_t z, slong len, slong prec)
{
    if (len == 0)
    {
        arb_poly_zero(res);
        return;
    }

    if (z->length < 2)
        len = 1;

    arb_poly_fit_length(res, len);
    _arb_hypgeom_bessel_j_series(res->coeffs, nu, z->coeffs, z->length, len, prec);
    _arb_poly_set_length(res, len);
    _arb_poly_normalise(res);
}
