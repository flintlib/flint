/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_dirichlet.h"

void
_acb_poly_zeta_cpx_series(acb_ptr z, const acb_t s, const acb_t a, int deflate, slong d, slong prec)
{
    ulong M, N;
    slong i, bound_prec;
    mag_t bound;
    arb_ptr vb;
    int is_real, const_is_real;

    if (d < 1)
        return;

    if (!acb_is_finite(s) || !acb_is_finite(a))
    {
        _acb_vec_indeterminate(z, d);
        return;
    }

    if (acb_is_one(s) && deflate && d == 1)
    {
        acb_digamma(z, a, prec);
        acb_neg(z, z);
        if (!acb_is_finite(z))  /* todo: in digamma */
            acb_indeterminate(z);
        return;
    }

    is_real = const_is_real = 0;

    if (acb_is_real(s) && acb_is_real(a))
    {
        if (arb_is_positive(acb_realref(a)))
        {
            is_real = const_is_real = 1;
        }
        else if (arb_is_int(acb_realref(a)) &&
             arb_is_int(acb_realref(s)) &&
             arb_is_nonpositive(acb_realref(s)))
        {
            const_is_real = 1;
        }
    }

    mag_init(bound);
    vb = _arb_vec_init(d);

    bound_prec = 40 + prec / 20;

    _acb_poly_zeta_em_choose_param(bound, &N, &M, s, a, FLINT_MIN(d, 2), prec, bound_prec);
    _acb_poly_zeta_em_bound(vb, s, a, N, M, d, bound_prec);

    _acb_poly_zeta_em_sum(z, s, a, deflate, N, M, d, prec);

    for (i = 0; i < d; i++)
    {
        arb_get_mag(bound, vb + i);
        arb_add_error_mag(acb_realref(z + i), bound);

        if (!is_real && !(i == 0 && const_is_real))
            arb_add_error_mag(acb_imagref(z + i), bound);
    }

    mag_clear(bound);
    _arb_vec_clear(vb, d);
}

void
_acb_poly_zeta_series(acb_ptr res, acb_srcptr h, slong hlen, const acb_t a, int deflate, slong len, slong prec)
{
    acb_ptr t, u;
    hlen = FLINT_MIN(hlen, len);

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);

    if (acb_is_one(a))
        acb_dirichlet_zeta_jet(t, h, deflate, len, prec);
    else
        _acb_poly_zeta_cpx_series(t, h, a, deflate, len, prec);

    /* compose with nonconstant part */
    acb_zero(u);
    _acb_vec_set(u + 1, h + 1, hlen - 1);
    _acb_poly_compose_series(res, t, len, u, hlen, len, prec);

    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
}

void
acb_poly_zeta_series(acb_poly_t res, const acb_poly_t f, const acb_t a, int deflate, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, n);

    if (f->length == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_poly_zeta_series(res->coeffs, t, 1, a, deflate, n, prec);
        acb_clear(t);
    }
    else
    {
        _acb_poly_zeta_series(res->coeffs, f->coeffs, f->length, a, deflate, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

