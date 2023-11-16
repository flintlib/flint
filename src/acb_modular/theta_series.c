/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_modular.h"

void
_acb_modular_theta_series(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4,
    acb_srcptr z, slong zlen, const acb_t tau, slong len, slong prec)
{
    acb_ptr t1, t2, t3, t4, t, v;

    zlen = FLINT_MIN(zlen, len);

    if (zlen <= 0)
        return;

    t = _acb_vec_init(4 * len);
    t1 = t;
    t2 = t1 + len;
    t3 = t2 + len;
    t4 = t3 + len;

    acb_modular_theta_jet(t1, t2, t3, t4, z, tau, len, prec);

    if (len == 1)
    {
        if (theta1 != NULL) acb_set(theta1, t1);
        if (theta2 != NULL) acb_set(theta2, t2);
        if (theta3 != NULL) acb_set(theta3, t3);
        if (theta4 != NULL) acb_set(theta4, t4);
    }
    else
    {
        v = _acb_vec_init(zlen);

        /* compose with nonconstant part */
        acb_zero(v);
        _acb_vec_set(v + 1, z + 1, zlen - 1);

        if (theta1 != NULL) _acb_poly_compose_series(theta1, t1, len, v, zlen, len, prec);
        if (theta2 != NULL) _acb_poly_compose_series(theta2, t2, len, v, zlen, len, prec);
        if (theta3 != NULL) _acb_poly_compose_series(theta3, t3, len, v, zlen, len, prec);
        if (theta4 != NULL) _acb_poly_compose_series(theta4, t4, len, v, zlen, len, prec);

        _acb_vec_clear(v, zlen);
    }

    _acb_vec_clear(t, 4 * len);
}

void
acb_modular_theta_series(acb_poly_t theta1, acb_poly_t theta2,
    acb_poly_t theta3, acb_poly_t theta4, const acb_poly_t z, const acb_t tau,
        slong len, slong prec)
{
    if (len == 0)
    {
        if (theta1 != NULL) acb_poly_zero(theta1);
        if (theta2 != NULL) acb_poly_zero(theta2);
        if (theta3 != NULL) acb_poly_zero(theta3);
        if (theta4 != NULL) acb_poly_zero(theta4);
        return;
    }

    if (z->length <= 1)
        len = 1;

    if (theta1 != NULL) acb_poly_fit_length(theta1, len);
    if (theta2 != NULL) acb_poly_fit_length(theta2, len);
    if (theta3 != NULL) acb_poly_fit_length(theta3, len);
    if (theta4 != NULL) acb_poly_fit_length(theta4, len);

    if (z->length == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_modular_theta_series(theta1 ? theta1->coeffs : NULL,
            theta2 ? theta2->coeffs : NULL, theta3 ? theta3->coeffs : NULL,
                theta4 ? theta4->coeffs : NULL, t, 1, tau, len, prec);
        acb_clear(t);
    }
    else
    {
        _acb_modular_theta_series(theta1 ? theta1->coeffs : NULL,
            theta2 ? theta2->coeffs : NULL, theta3 ? theta3->coeffs : NULL,
                theta4 ? theta4->coeffs : NULL, z->coeffs, z->length, tau, len, prec);
    }

    if (theta1 != NULL) _acb_poly_set_length(theta1, len);
    if (theta2 != NULL) _acb_poly_set_length(theta2, len);
    if (theta3 != NULL) _acb_poly_set_length(theta3, len);
    if (theta4 != NULL) _acb_poly_set_length(theta4, len);

    if (theta1 != NULL) _acb_poly_normalise(theta1);
    if (theta2 != NULL) _acb_poly_normalise(theta2);
    if (theta3 != NULL) _acb_poly_normalise(theta3);
    if (theta4 != NULL) _acb_poly_normalise(theta4);
}

