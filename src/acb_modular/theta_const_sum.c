/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_modular.h"

void
acb_modular_theta_const_sum(acb_t theta2, acb_t theta3, acb_t theta4,
    const acb_t q, slong prec)
{
    mag_t qmag, err;
    double log2q_approx;
    int is_real, is_real_or_imag;
    slong N;

    mag_init(qmag);
    mag_init(err);

    acb_get_mag(qmag, q);
    log2q_approx = mag_get_d_log2_approx(qmag);

    is_real = arb_is_zero(acb_imagref(q));
    is_real_or_imag = is_real || arb_is_zero(acb_realref(q));

    if (log2q_approx >= 0.0)
    {
        N = 1;
        mag_inf(err);
    }
    else
    {
        N = 0;

        while (0.05 * N * N < prec)
        {
            if (log2q_approx * ((N+2)*(N+2)/4) < -prec - 2)
                break;
            N++;
        }
        N = (N+2)*(N+2)/4;

        mag_geom_series(err, qmag, N);
        mag_mul_2exp_si(err, err, 1); /* each term is taken twice */

        if (mag_is_inf(err))
            N = 1;
    }

    if (N < 1800)
        acb_modular_theta_const_sum_basecase(theta2, theta3, theta4, q, N, prec);
    else
        acb_modular_theta_const_sum_rs(theta2, theta3, theta4, q, N, prec);

    if (is_real_or_imag)
        arb_add_error_mag(acb_realref(theta2), err);
    else
        acb_add_error_mag(theta2, err);

    if (is_real)
    {
        arb_add_error_mag(acb_realref(theta3), err);
        arb_add_error_mag(acb_realref(theta4), err);
    }
    else
    {
        acb_add_error_mag(theta3, err);
        acb_add_error_mag(theta4, err);
    }

    mag_clear(qmag);
    mag_clear(err);
}

