/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void _arb_poly_lambertw_series(arb_ptr res,
    arb_srcptr z, slong zlen, int flags, slong len, slong prec)
{
    arb_ptr w, ew, t, u;
    arb_t ew0;

    zlen = FLINT_MIN(zlen, len);

    if (zlen == 1)
    {
        arb_lambertw(res, z, flags, prec);
        _arb_vec_zero(res + 1, len - 1);
        return;
    }

    w = _arb_vec_init(len);
    ew = _arb_vec_init(len);
    t = _arb_vec_init(len);
    u = _arb_vec_init(len);
    arb_init(ew0);

    arb_lambertw(w, z, flags, prec);

    if (arb_contains_zero(w))
        arb_exp(ew0, w, prec);
    else
        arb_div(ew0, z, w, prec);

    arb_add(t, ew0, z, prec);
    arb_div(w + 1, z + 1, t, prec);

    NEWTON_INIT(2, len)
    NEWTON_LOOP(m, n)

    /* _arb_poly_exp_series(ew, w, m, n, prec); */
    arb_zero(t);
    _arb_vec_set(t + 1, w + 1, m - 1);
    _arb_poly_exp_series(ew, t, m, n, prec);
    _arb_vec_scalar_mul(ew, ew, n, ew0, prec);

    _arb_poly_mullow(t, ew, n, w, m, n, prec);
    _arb_poly_sub(u, t, n, z, FLINT_MIN(zlen, n), prec);
    _arb_vec_add(t, t, ew, n, prec);
    _arb_poly_div_series(ew, u, n, t, n, n, prec);
    _arb_vec_neg(w + m, ew + m, n - m);

    NEWTON_END_LOOP
    NEWTON_END

    _arb_vec_set(res, w, len);

    _arb_vec_clear(w, len);
    _arb_vec_clear(ew, len);
    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);
    arb_clear(ew0);
}

void
arb_poly_lambertw_series(arb_poly_t res,
    const arb_poly_t z, int flags, slong len, slong prec)
{
    if (len == 0 || (flags == 0 && z->length == 0))
    {
        arb_poly_zero(res);
        return;
    }

    if (z->length == 0)
    {
        arb_poly_fit_length(res, len);
        _arb_vec_indeterminate(res->coeffs, len);
        _arb_poly_set_length(res, len);
        return;
    }

    arb_poly_fit_length(res, len);
    _arb_poly_lambertw_series(res->coeffs, z->coeffs, z->length, flags, len, prec);
    _arb_poly_set_length(res, len);
    _arb_poly_normalise(res);
}

