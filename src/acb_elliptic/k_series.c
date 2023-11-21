/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_elliptic.h"

void
_acb_elliptic_k_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
{
    acb_ptr t, u;

    zlen = FLINT_MIN(zlen, len);

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);

    acb_elliptic_k_jet(t, z, len, prec);

    /* compose with nonconstant part */
    acb_zero(u);
    _acb_vec_set(u + 1, z + 1, zlen - 1);
    _acb_poly_compose_series(res, t, len, u, zlen, len, prec);

    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
}

void
acb_elliptic_k_series(acb_poly_t res, const acb_poly_t z, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, n);

    if (z->length == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_elliptic_k_series(res->coeffs, t, 1, n, prec);
        acb_clear(t);
    }
    else
    {
        _acb_elliptic_k_series(res->coeffs, z->coeffs, z->length, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}
