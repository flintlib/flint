/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

void _acb_poly_lambertw_series(acb_ptr res,
    acb_srcptr z, slong zlen, const fmpz_t k, int flags, slong len, slong prec)
{
    acb_ptr w, ew, t, u;
    acb_t ew0;

    zlen = FLINT_MIN(zlen, len);

    if (zlen == 1)
    {
        acb_lambertw(res, z, k, flags, prec);
        _acb_vec_zero(res + 1, len - 1);
        return;
    }

    w = _acb_vec_init(len);
    ew = _acb_vec_init(len);
    t = _acb_vec_init(len);
    u = _acb_vec_init(len);
    acb_init(ew0);

    acb_lambertw(w, z, k, flags, prec);

    if (acb_contains_zero(w))
        acb_exp(ew0, w, prec);
    else
        acb_div(ew0, z, w, prec);

    acb_add(t, ew0, z, prec);
    acb_div(w + 1, z + 1, t, prec);

    NEWTON_INIT(2, len)
    NEWTON_LOOP(m, n)

    /* _acb_poly_exp_series(ew, w, m, n, prec); */
    acb_zero(t);
    _acb_vec_set(t + 1, w + 1, m - 1);
    _acb_poly_exp_series(ew, t, m, n, prec);
    _acb_vec_scalar_mul(ew, ew, n, ew0, prec);

    _acb_poly_mullow(t, ew, n, w, m, n, prec);
    _acb_poly_sub(u, t, n, z, FLINT_MIN(zlen, n), prec);
    _acb_vec_add(t, t, ew, n, prec);
    _acb_poly_div_series(ew, u, n, t, n, n, prec);
    _acb_vec_neg(w + m, ew + m, n - m);

    NEWTON_END_LOOP
    NEWTON_END

    _acb_vec_set(res, w, len);

    _acb_vec_clear(w, len);
    _acb_vec_clear(ew, len);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
    acb_clear(ew0);
}

void
acb_poly_lambertw_series(acb_poly_t res,
    const acb_poly_t z, const fmpz_t k, int flags, slong len, slong prec)
{
    if (len == 0 || (fmpz_is_zero(k) && z->length == 0))
    {
        acb_poly_zero(res);
        return;
    }

    if (z->length == 0)
    {
        acb_poly_fit_length(res, len);
        _acb_vec_indeterminate(res->coeffs, len);
        _acb_poly_set_length(res, len);
        return;
    }

    acb_poly_fit_length(res, len);
    _acb_poly_lambertw_series(res->coeffs, z->coeffs, z->length, k, flags, len, prec);
    _acb_poly_set_length(res, len);
    _acb_poly_normalise(res);
}

