/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_hypgeom.h"

void
_arb_hypgeom_beta_lower_series(arb_ptr res,
    const arb_t a, const arb_t b, arb_srcptr z, slong zlen, int regularized,
    slong len, slong prec)
{
    arb_ptr t, u, v;
    arb_t c, d, e;

    zlen = FLINT_MIN(zlen, len);

    if (zlen == 1)
    {
        arb_hypgeom_beta_lower(res, a, b, z, regularized, prec);
        _arb_vec_zero(res + 1, len - 1);
        return;
    }

    t = _arb_vec_init(len);
    u = _arb_vec_init(len);
    v = _arb_vec_init(zlen - 1);
    arb_init(c);
    arb_init(d);
    arb_init(e);

    arb_hypgeom_beta_lower(d, a, b, z, regularized, prec);

    if (regularized)
    {
        /* todo: except in special cases, we already compute a bunch of
           gamma functions in beta_lower, so we could avoid recomputing them */
        arb_add(e, a, b, prec);
        arb_gamma(e, e, prec);
        arb_rgamma(c, a, prec);
        arb_mul(e, e, c, prec);
        arb_rgamma(c, b, prec);
        arb_mul(e, e, c, prec);
    }

    /* u = (1-z)^(b-1) */
    _arb_vec_neg(t, z, zlen);
    arb_add_ui(t, t, 1, prec);
    arb_sub_ui(c, b, 1, prec);
    _arb_poly_pow_arb_series(u, t, FLINT_MIN(zlen, len - 1), c, len - 1, prec);

    /* t = z^(a-1) */
    arb_sub_ui(c, a, 1, prec);
    _arb_poly_pow_arb_series(t, z, FLINT_MIN(zlen, len - 1), c, len - 1, prec);

    /* v = z' */
    _arb_poly_derivative(v, z, zlen, prec);

    _arb_poly_mullow(res, t, len - 1, u, len - 1, len - 1, prec);
    _arb_poly_mullow(t, res, len - 1, v, zlen - 1, len - 1, prec);

    _arb_poly_integral(res, t, len, prec);

    if (regularized)
    {
        _arb_vec_scalar_mul(res, res, len, e, prec);
    }

    arb_set(res, d);

    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);
    _arb_vec_clear(v, zlen - 1);
    arb_clear(c);
    arb_clear(d);
    arb_clear(e);
}

void arb_hypgeom_beta_lower_series(arb_poly_t res,
    const arb_t a, const arb_t b, const arb_poly_t z, int regularized,
    slong len, slong prec)
{
    if (len == 0)
    {
        arb_poly_zero(res);
        return;
    }

    arb_poly_fit_length(res, len);

    if (z->length == 0)
    {
        arb_t t;
        arb_init(t);
        _arb_hypgeom_beta_lower_series(res->coeffs, a, b, t, 1,
            regularized, len, prec);
        arb_clear(t);
    }
    else
    {
        _arb_hypgeom_beta_lower_series(res->coeffs, a, b, z->coeffs,
            z->length, regularized, len, prec);
    }

    _arb_poly_set_length(res, len);
    _arb_poly_normalise(res);
}

