/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"
#include "acb_modular.h"

void
acb_elliptic_p_prime(acb_t r, const acb_t z, const acb_t tau, slong prec)
{
    acb_struct tz[4];
    acb_t t1, t2, t3;
    int i, real;

    real = acb_is_real(z) && arb_is_int_2exp_si(acb_realref(tau), -1) &&
                arb_is_positive(acb_imagref(tau));

    acb_init(t1);
    acb_init(t2);
    acb_init(t3);

    for (i = 0; i < 4; i++)
        acb_init(tz + i);

    acb_modular_theta(tz + 0, tz + 1, tz + 2, tz + 3, z, tau, prec);

    /* (-2*pi*eta^2/theta1)^3*theta2*theta3*theta4 */
    acb_const_pi(t2, prec);
    acb_mul_2exp_si(t2, t2, 1);
    acb_neg(t2, t2);
    acb_modular_eta(t3, tau, prec);
    acb_mul(t1, t3, t3, prec);
    acb_mul(t3, t1, t2, prec);
    acb_div(t1, t3, tz + 0, prec);
    acb_mul(t2, t1, t1, prec);
    acb_mul(t3, t1, t2, prec);
    acb_mul(t1, tz + 1, tz + 2, prec);
    acb_mul(t2, t1, tz + 3, prec);
    acb_mul(r, t3, t2, prec);

    if (real)
        arb_zero(acb_imagref(r));

    acb_clear(t1);
    acb_clear(t2);
    acb_clear(t3);

    for (i = 0; i < 4; i++)
        acb_clear(tz + i);
}

