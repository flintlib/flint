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

/* series of c^(d+x) */
static inline void
_arb_poly_pow_cpx(arb_ptr res, const arb_t c, const arb_t d, slong trunc, slong prec)
{
    slong i;
    arb_t logc;

    arb_init(logc);
    arb_log(logc, c, prec);
    arb_mul(res + 0, logc, d, prec);
    arb_exp(res + 0, res + 0, prec);

    for (i = 1; i < trunc; i++)
    {
        arb_mul(res + i, res + i - 1, logc, prec);
        arb_div_ui(res + i, res + i, i, prec);
    }

    arb_clear(logc);
}

void
_arb_poly_zeta_series(arb_ptr res, arb_srcptr h, slong hlen, const arb_t a, int deflate, slong len, slong prec)
{
    slong i;
    acb_t cs, ca;
    acb_ptr z;
    arb_ptr t, u;

    if (arb_contains_nonpositive(a))
    {
        _arb_vec_indeterminate(res, len);
        return;
    }

    hlen = FLINT_MIN(hlen, len);

    z = _acb_vec_init(len);
    t = _arb_vec_init(len);
    u = _arb_vec_init(len);
    acb_init(cs);
    acb_init(ca);

    /* use reflection formula */
    if (arf_sgn(arb_midref(h)) < 0 && arb_is_one(a))
    {
        /* zeta(s) = (2*pi)**s * sin(pi*s/2) / pi * gamma(1-s) * zeta(1-s) */
        arb_t pi;
        arb_ptr f, s1, s2, s3, s4;

        arb_init(pi);
        f = _arb_vec_init(2);
        s1 = _arb_vec_init(len);
        s2 = _arb_vec_init(len);
        s3 = _arb_vec_init(len);
        s4 = _arb_vec_init(len);

        arb_const_pi(pi, prec);

        /* s1 = (2*pi)**s */
        arb_mul_2exp_si(pi, pi, 1);
        _arb_poly_pow_cpx(s1, pi, h, len, prec);
        arb_mul_2exp_si(pi, pi, -1);

        /* s2 = sin(pi*s/2) / pi */
        arb_set(f, h);
        arb_one(f + 1);
        arb_mul_2exp_si(f, f, -1);
        arb_mul_2exp_si(f + 1, f + 1, -1);
        _arb_poly_sin_pi_series(s2, f, 2, len, prec);
        _arb_vec_scalar_div(s2, s2, len, pi, prec);

        /* s3 = gamma(1-s) */
        arb_sub_ui(f, h, 1, prec);
        arb_neg(f, f);
        arb_set_si(f + 1, -1);
        _arb_poly_gamma_series(s3, f, 2, len, prec);

        /* s4 = zeta(1-s) */
        arb_sub_ui(f, h, 1, prec);
        arb_neg(f, f);
        acb_set_arb(cs, f);
        acb_one(ca);
        _acb_poly_zeta_cpx_series(z, cs, ca, 0, len, prec);
        for (i = 0; i < len; i++)
            arb_set(s4 + i, acb_realref(z + i));
        for (i = 1; i < len; i += 2)
            arb_neg(s4 + i, s4 + i);

        _arb_poly_mullow(u, s1, len, s2, len, len, prec);
        _arb_poly_mullow(s1, s3, len, s4, len, len, prec);
        _arb_poly_mullow(t, u, len, s1, len, len, prec);

        /* add 1/(1-(s+t)) = 1/(1-s) + t/(1-s)^2 + ... */
        if (deflate)
        {
            arb_sub_ui(u, h, 1, prec);
            arb_neg(u, u);
            arb_inv(u, u, prec);
            for (i = 1; i < len; i++)
                arb_mul(u + i, u + i - 1, u, prec);
            _arb_vec_add(t, t, u, len, prec);
        }

        arb_clear(pi);
        _arb_vec_clear(f, 2);
        _arb_vec_clear(s1, len);
        _arb_vec_clear(s2, len);
        _arb_vec_clear(s3, len);
        _arb_vec_clear(s4, len);
    }
    else
    {
        acb_set_arb(cs, h);
        acb_set_arb(ca, a);
        _acb_poly_zeta_cpx_series(z, cs, ca, deflate, len, prec);
        for (i = 0; i < len; i++)
            arb_set(t + i, acb_realref(z + i));
    }

    /* compose with nonconstant part */
    arb_zero(u);
    _arb_vec_set(u + 1, h + 1, hlen - 1);
    _arb_poly_compose_series(res, t, len, u, hlen, len, prec);

    _acb_vec_clear(z, len);
    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);
    acb_clear(cs);
    acb_clear(ca);
}

void
arb_poly_zeta_series(arb_poly_t res, const arb_poly_t f, const arb_t a, int deflate, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    arb_poly_fit_length(res, n);

    if (f->length == 0)
    {
        arb_t t;
        arb_init(t);
        _arb_poly_zeta_series(res->coeffs, t, 1, a, deflate, n, prec);
        arb_clear(t);
    }
    else
    {
        _arb_poly_zeta_series(res->coeffs, f->coeffs, f->length, a, deflate, n, prec);
    }

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

