/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_hypgeom_beta_lower_series(acb_ptr res,
    const acb_t a, const acb_t b, acb_srcptr z, slong zlen, int regularized,
    slong len, slong prec)
{
    acb_ptr t, u, v;
    acb_t c, d, e;

    zlen = FLINT_MIN(zlen, len);

    if (zlen == 1)
    {
        acb_hypgeom_beta_lower(res, a, b, z, regularized, prec);
        _acb_vec_zero(res + 1, len - 1);
        return;
    }

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);
    v = _acb_vec_init(zlen - 1);
    acb_init(c);
    acb_init(d);
    acb_init(e);

    acb_hypgeom_beta_lower(d, a, b, z, regularized, prec);

    if (regularized)
    {
        /* todo: except in special cases, we already compute a bunch of
           gamma functions in beta_lower, so we could avoid recomputing them */
        acb_add(e, a, b, prec);
        acb_gamma(e, e, prec);
        acb_rgamma(c, a, prec);
        acb_mul(e, e, c, prec);
        acb_rgamma(c, b, prec);
        acb_mul(e, e, c, prec);
    }

    /* u = (1-z)^(b-1) */
    _acb_vec_neg(t, z, zlen);
    acb_add_ui(t, t, 1, prec);
    acb_sub_ui(c, b, 1, prec);
    _acb_poly_pow_acb_series(u, t, FLINT_MIN(zlen, len - 1), c, len - 1, prec);

    /* t = z^(a-1) */
    acb_sub_ui(c, a, 1, prec);
    _acb_poly_pow_acb_series(t, z, FLINT_MIN(zlen, len - 1), c, len - 1, prec);

    /* v = z' */
    _acb_poly_derivative(v, z, zlen, prec);

    _acb_poly_mullow(res, t, len - 1, u, len - 1, len - 1, prec);
    _acb_poly_mullow(t, res, len - 1, v, zlen - 1, len - 1, prec);

    _acb_poly_integral(res, t, len, prec);

    if (regularized)
    {
        _acb_vec_scalar_mul(res, res, len, e, prec);
    }

    acb_set(res, d);

    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
    _acb_vec_clear(v, zlen - 1);
    acb_clear(c);
    acb_clear(d);
    acb_clear(e);
}

void acb_hypgeom_beta_lower_series(acb_poly_t res,
    const acb_t a, const acb_t b, const acb_poly_t z, int regularized,
    slong len, slong prec)
{
    if (len == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, len);

    if (z->length == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_hypgeom_beta_lower_series(res->coeffs, a, b, t, 1,
            regularized, len, prec);
        acb_clear(t);
    }
    else
    {
        _acb_hypgeom_beta_lower_series(res->coeffs, a, b, z->coeffs,
            z->length, regularized, len, prec);
    }

    _acb_poly_set_length(res, len);
    _acb_poly_normalise(res);
}

