/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/*
f(s)  = (f(s+hi) + f(s-hi))/2     + err0
f'(s) = (f(s+hi) - f(s-hi))/(2hi) + err1

|err0| <= |f(s +/- R)| (1/R)^2 h^2 / (1 - h/R)
|err1| <= |f(s +/- R)| (1/R)^3 h^2 / (1 - h/R)

assuming h/R < 1
*/
static void
acb_dirichlet_zeta_jet_rs_mid(acb_ptr res, const acb_t s, slong prec)
{
    acb_t t, u;
    arb_t hh;
    mag_t h, R, err0, err1, tmp;
    slong hexp;

    acb_init(t);
    acb_init(u);
    arb_init(hh);
    mag_init(h);
    mag_init(R);
    mag_init(err0);
    mag_init(err1);
    mag_init(tmp);

    hexp = -(prec / 2) - 10;
    mag_set_ui_2exp_si(h, 1, hexp);
    mag_set_ui_2exp_si(R, 1, -5);

    acb_set(t, s);
    mag_add(arb_radref(acb_realref(t)), arb_radref(acb_realref(t)), R);
    mag_add(arb_radref(acb_imagref(t)), arb_radref(acb_imagref(t)), R);
    /* tmp = |f(s +/- R)| */
    acb_dirichlet_zeta_bound(tmp, t);

    /* err0 = (1/R)^2 h^2 / (1 - h/R) */
    mag_div(err0, h, R);
    mag_geom_series(err0, err0, 2);
    /* err0 *= tmp */
    mag_mul(err0, err0, tmp);
    /* err1 = err0 / R */
    mag_div(err1, err0, R);

    /* hh = h as an arb_t */
    arb_one(hh);
    arb_mul_2exp_si(hh, hh, hexp);

    acb_set(t, s);
    acb_set(u, s);
    arb_add(acb_imagref(t), acb_imagref(t), hh, 10 * prec);
    arb_sub(acb_imagref(u), acb_imagref(u), hh, 10 * prec);

    /* zeta(mid(s)+ih) */
    acb_dirichlet_zeta_rs(res, t, 0, 1.5 * prec + 10);
    /* zeta(mid(s)-ih) */
    acb_dirichlet_zeta_rs(res + 1, u, 0, 1.5 * prec + 10);

    acb_sub(t, res, res + 1, prec);
    acb_add(res, res, res + 1, prec);
    acb_swap(res + 1, t);

    acb_mul_2exp_si(res, res, -1);
    acb_mul_2exp_si(res + 1, res + 1, -1);
    acb_mul_2exp_si(res + 1, res + 1, -hexp);
    acb_div_onei(res + 1, res + 1);

    acb_add_error_mag(res, err0);
    acb_add_error_mag(res + 1, err1);

    acb_clear(t);
    acb_clear(u);
    arb_clear(hh);
    mag_clear(h);
    mag_clear(R);
    mag_clear(err0);
    mag_clear(err1);
    mag_clear(tmp);
}

void
acb_dirichlet_zeta_jet_rs(acb_ptr res, const acb_t s, slong len, slong prec)
{
    if (len > 2)
    {
        flint_throw(FLINT_ERROR, "acb_dirichlet_zeta_jet_rs: len > 2 not implemented\n");
    }

    if (len <= 0)
        return;

    if (len == 1)
    {
        acb_dirichlet_zeta_rs(res, s, 0, prec);
    }
    else if (acb_is_exact(s))
    {
        acb_dirichlet_zeta_jet_rs_mid(res, s, prec);
    }
    else
    {
        acb_t t;
        mag_t r, err0, err1, der1, der2, M;

        /*
        slong acc;

        acc = acb_rel_accuracy_bits(s);
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + 20);
        */

        /*
        assume s = m +/- r

        f(s)  = f(m) + err0
        f'(s) = f'(m) + err1

        |err1| <= |f''(s)| r
        |err0| <= min(|f'(s)| r, |f'(m)| r + 0.5 |f''(s)| r^2)
                = r min(|f'(s)|, |f'(m)| + 0.5 |f''(s)| r)
        */

        acb_init(t);
        mag_init(r);
        mag_init(err0);
        mag_init(err1);
        mag_init(der1);
        mag_init(der2);
        mag_init(M);

        /* r = rad(s) */
        mag_hypot(r, arb_radref(acb_realref(s)), arb_radref(acb_imagref(s)));

        /* Bound zeta'(s), zeta''(s) */
        acb_dirichlet_zeta_deriv_bound(der1, der2, s);

        /* f(m), f'(m) */
        acb_get_mid(t, s);
        acb_dirichlet_zeta_jet_rs_mid(res, t, prec);

        /* err1 = |f''(s)| r */
        mag_mul(err1, der2, r);

        /* err0 = |f'(m)| + 0.5 |f''(s)| r */
        acb_get_mag(M, res + 1);
        mag_mul_2exp_si(err0, err1, -1);
        mag_add(err0, err0, M);
        /* err0 = min(err0, |f'(s)| */
        mag_min(err0, err0, der1);
        /* err0 = err0 * r */
        mag_mul(err0, err0, r);

        acb_add_error_mag(res, err0);
        acb_add_error_mag(res + 1, err1);

        acb_clear(t);
        mag_clear(r);
        mag_clear(err0);
        mag_clear(err1);
        mag_clear(der1);
        mag_clear(der2);
        mag_clear(M);
    }
}

