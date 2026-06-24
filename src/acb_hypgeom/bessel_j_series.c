/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"
#include "gr_poly.h"

void
_acb_hypgeom_bessel_j_series(acb_ptr res, const acb_t nu,
    acb_srcptr z, slong zlen, slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    if (_gr_poly_bessel_j_series(res, nu, z, zlen, len, ctx) != GR_SUCCESS)
        _acb_vec_indeterminate(res, len);
}

void
acb_hypgeom_bessel_j_series(acb_poly_t res, const acb_t nu,
    const acb_poly_t z, slong len, slong prec)
{
    if (len == 0)
    {
        acb_poly_zero(res);
        return;
    }

    if (z->length < 2)
        len = 1;

    acb_poly_fit_length(res, len);
    _acb_hypgeom_bessel_j_series(res->coeffs, nu, z->coeffs, z->length, len, prec);
    _acb_poly_set_length(res, len);
    _acb_poly_normalise(res);
}
