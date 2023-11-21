/*
    Copyright (C) 2014 Fredrik Johansson

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
acb_elliptic_p(acb_t r, const acb_t z, const acb_t tau, slong prec)
{
    acb_struct t0[4], tz[4];
    acb_t t;
    int i, real;

    real = acb_is_real(z) && arb_is_int_2exp_si(acb_realref(tau), -1) &&
                arb_is_positive(acb_imagref(tau));

    acb_init(t);

    for (i = 0; i < 4; i++)
    {
        acb_init(t0 + i);
        acb_init(tz + i);
    }

    acb_modular_theta(tz + 0, tz + 1, tz + 2, tz + 3, z, tau, prec);
    acb_zero(t);
    acb_modular_theta(t0 + 0, t0 + 1, t0 + 2, t0 + 3, t, tau, prec);

    acb_mul(t, t0 + 1, t0 + 2, prec);
    acb_mul(t, t, tz + 3, prec);
    acb_div(t, t, tz + 0, prec);
    acb_mul(t, t, t, prec);

    acb_pow_ui(t0 + 1, t0 + 1, 4, prec);
    acb_pow_ui(t0 + 2, t0 + 2, 4, prec);
    acb_add(r, t0 + 1, t0 + 2, prec);
    acb_div_ui(r, r, 3, prec);

    acb_sub(r, t, r, prec);

    acb_const_pi(t, prec);
    acb_mul(t, t, t, prec);
    acb_mul(r, r, t, prec);

    if (real)
        arb_zero(acb_imagref(r));

    acb_clear(t);

    for (i = 0; i < 4; i++)
    {
        acb_clear(t0 + i);
        acb_clear(tz + i);
    }
}

