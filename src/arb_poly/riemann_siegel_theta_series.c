/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

void
_arb_poly_riemann_siegel_theta_series(arb_ptr res,
    arb_srcptr h, slong hlen, slong len, slong prec)
{
    acb_ptr s;
    arb_t u;
    slong i;

    hlen = FLINT_MIN(hlen, len);

    s = _acb_vec_init(len);

    arb_init(u);

    /* s = 1/4 + (1/2) i h */
    for (i = 0; i < hlen; i++)
        arb_mul_2exp_si(acb_imagref(s + i), h + i, -1);

    arb_one(u);
    arb_mul_2exp_si(u, u, -2);
    arb_add(acb_realref(s), acb_realref(s), u, prec);

    /* log gamma */
    _acb_poly_lgamma_series(s, s, hlen, len, prec);

    /* imaginary part */
    for (i = 0; i < len; i++)
        arb_set(res + i, acb_imagref(s + i));

    /* subtract log(pi)/2 * h */
    arb_const_pi(u, prec);
    arb_log(u, u, prec);
    arb_mul_2exp_si(u, u, -1);
    arb_neg(u, u);
    _arb_vec_scalar_addmul(res, h, hlen, u, prec);

    _acb_vec_clear(s, len);
    arb_clear(u);
}

void
arb_poly_riemann_siegel_theta_series(arb_poly_t res,
    const arb_poly_t f, slong n, slong prec)
{
    if (n == 0 || f->length == 0)
    {
        arb_poly_zero(res);
        return;
    }

    if (res == f)
    {
        arb_poly_t tmp;
        arb_poly_init2(tmp, n);
        _arb_poly_riemann_siegel_theta_series(tmp->coeffs,
            f->coeffs, f->length, n, prec);
        arb_poly_swap(res, tmp);
        arb_poly_clear(tmp);
    }
    else
    {
        arb_poly_fit_length(res, n);
        _arb_poly_riemann_siegel_theta_series(res->coeffs,
            f->coeffs, f->length, n, prec);
    }

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

