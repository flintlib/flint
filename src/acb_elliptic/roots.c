/*
    Copyright (C) 2017 Fredrik Johansson

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
acb_elliptic_roots(acb_t e1, acb_t e2, acb_t e3, const acb_t tau, slong prec)
{
    acb_t t1, t2, t3, t4;
    int e1r, e23r;

    if (!arb_is_positive(acb_imagref(tau)) || !arb_is_finite(acb_realref(tau)))
    {
        acb_indeterminate(e1);
        acb_indeterminate(e2);
        acb_indeterminate(e3);
        return;
    }

    acb_init(t1);
    acb_init(t2);
    acb_init(t3);
    acb_init(t4);

    e1r = e23r = 0;

    if (arb_is_int(acb_realref(tau)))
        e1r = e23r = 1;
    else if (arb_is_int_2exp_si(acb_realref(tau), -1))
        e1r = 1;

    acb_modular_theta(t1, t2, t3, t4, t1, tau, prec);

    acb_pow_ui(t2, t2, 4, prec);
    acb_pow_ui(t4, t4, 4, prec);

    acb_sub(e2, t2, t4, prec);
    acb_mul_2exp_si(t3, t4, 1);
    acb_add(e1, t2, t3, prec);
    acb_mul_2exp_si(t3, t2, 1);
    acb_add(e3, t3, t4, prec);

    acb_const_pi(t3, prec);
    acb_mul(t3, t3, t3, prec);
    acb_div_ui(t3, t3, 3, prec);

    acb_mul(e1, e1, t3, prec);
    acb_mul(e2, e2, t3, prec);
    acb_mul(e3, e3, t3, prec);
    acb_neg(e3, e3);

    if (e1r)
        arb_zero(acb_imagref(e1));

    if (e23r)
    {
        arb_zero(acb_imagref(e2));
        arb_zero(acb_imagref(e3));
    }

    acb_clear(t1);
    acb_clear(t2);
    acb_clear(t3);
    acb_clear(t4);
}

