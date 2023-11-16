/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_si_asymp(acb_t res, const acb_t z, slong prec)
{
    acb_t t, u, w, v, one;

    acb_init(t);
    acb_init(u);
    acb_init(w);
    acb_init(v);
    acb_init(one);

    acb_one(one);
    acb_mul_onei(w, z);

    /* u = U(1,1,iz) */
    acb_hypgeom_u_asymp(u, one, one, w, -1, prec);
    /* v = e^(-iz) */
    acb_neg(v, w);
    acb_exp(v, v, prec);
    acb_mul(t, u, v, prec);

    if (acb_is_real(z))
    {
        arb_div(acb_realref(t), acb_realref(t), acb_realref(z), prec);
        arb_zero(acb_imagref(t));
        acb_neg(t, t);
    }
    else
    {
        /* u = U(1,1,-iz) */
        acb_neg(w, w);
        acb_hypgeom_u_asymp(u, one, one, w, -1, prec);
        acb_inv(v, v, prec);
        acb_addmul(t, u, v, prec);

        acb_div(t, t, z, prec);
        acb_mul_2exp_si(t, t, -1);
        acb_neg(t, t);
    }

    if (arb_is_zero(acb_realref(z)))
    {
        /* the function is imaginary */
        arb_zero(acb_realref(t));
    }
    else if (arb_is_positive(acb_realref(z)))
    {
        acb_const_pi(u, prec);
        acb_mul_2exp_si(u, u, -1);
        arb_add(acb_realref(t), acb_realref(t), acb_realref(u), prec);
    }
    else if (arb_is_negative(acb_realref(z)))
    {
        acb_const_pi(u, prec);
        acb_mul_2exp_si(u, u, -1);
        arb_sub(acb_realref(t), acb_realref(t), acb_realref(u), prec);
    }
    else
    {
        /* add [-pi,pi]/2 */
        acb_const_pi(u, prec);
        acb_mul_2exp_si(u, u, -1);
        arb_add_error(acb_imagref(t), acb_realref(u));
    }

    acb_swap(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(w);
    acb_clear(v);
    acb_clear(one);
}

void
acb_hypgeom_si_1f2(acb_t res, const acb_t z, slong prec)
{
    acb_t a, t;
    acb_struct b[3];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(b + 2);
    acb_init(t);

    acb_one(a);
    acb_mul_2exp_si(a, a, -1);
    acb_set_ui(b, 3);
    acb_mul_2exp_si(b, b, -1);
    acb_set(b + 1, b);
    acb_one(b + 2);

    acb_mul(t, z, z, prec);
    acb_mul_2exp_si(t, t, -2);
    acb_neg(t, t);
    acb_hypgeom_pfq_direct(t, a, 1, b, 3, t, -1, prec);
    acb_mul(t, t, z, prec);

    acb_swap(res, t);

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(b + 2);
    acb_clear(t);
}

void
acb_hypgeom_si(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_real(z) && arb_is_finite(acb_realref(z)))
    {
        arb_hypgeom_si(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
        return;
    }

    if (acb_hypgeom_u_use_asymp(z, prec))
        acb_hypgeom_si_asymp(res, z, prec);
    else
        acb_hypgeom_si_1f2(res, z, prec);
}
