/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_elliptic.h"
#include "acb_modular.h"
#include "acb_poly.h"

void
acb_elliptic_p_jet(acb_ptr r, const acb_t z, const acb_t tau, slong len, slong prec)
{
    acb_t t01, t02, t03, t04;
    acb_ptr tz1, tz2, tz3, tz4;
    acb_t t;
    int real;
    slong k;

    if (len < 1)
        return;

    if (len == 1)
    {
        acb_elliptic_p(r, z, tau, prec);
        return;
    }

    real = acb_is_real(z) && arb_is_int_2exp_si(acb_realref(tau), -1) &&
                arb_is_positive(acb_imagref(tau));

    acb_init(t);

    acb_init(t01);
    acb_init(t02);
    acb_init(t03);
    acb_init(t04);

    tz1 = _acb_vec_init(len);
    tz2 = _acb_vec_init(len);
    tz3 = _acb_vec_init(len);
    tz4 = _acb_vec_init(len);

    acb_modular_theta_jet(tz1, tz2, tz3, tz4, z, tau, len, prec);

    /* [theta_4(z) / theta_1(z)]^2 */
    _acb_poly_div_series(tz2, tz4, len, tz1, len, len, prec);
    _acb_poly_mullow(tz1, tz2, len, tz2, len, len, prec);

    acb_zero(t);
    acb_modular_theta(t01, t02, t03, t04, t, tau, prec);

    /* [theta_2(0) * theta_3(0)] ^2 */
    acb_mul(t, t02, t03, prec);
    acb_mul(t, t, t, prec);
    _acb_vec_scalar_mul(tz1, tz1, len, t, prec);

    /* - [theta_2(0)^4 + theta_3(0)^4] / 3 */
    acb_pow_ui(t02, t02, 4, prec);
    acb_pow_ui(t03, t03, 4, prec);
    acb_add(t, t02, t03, prec);
    acb_div_ui(t, t, 3, prec);
    acb_sub(tz1, tz1, t, prec);

    /* times pi^2 */
    acb_const_pi(t, prec);
    acb_mul(t, t, t, prec);
    _acb_vec_scalar_mul(r, tz1, len, t, prec);

    if (real)
    {
        for (k = 0; k < len; k++)
            arb_zero(acb_imagref(r + k));
    }

    acb_clear(t);

    acb_clear(t01);
    acb_clear(t02);
    acb_clear(t03);
    acb_clear(t04);

    _acb_vec_clear(tz1, len);
    _acb_vec_clear(tz2, len);
    _acb_vec_clear(tz3, len);
    _acb_vec_clear(tz4, len);
}

