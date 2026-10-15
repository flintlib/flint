/*
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2026 Ilia Gradina

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_hypgeom.h"

/* K_nu(z) = K_{nu+2}(z) - 2 (nu+1) K_{nu+1}(z) / z */
static void
bessel_k_integration_recurrence(arb_t res, const arb_t nu, const arb_t z,
    int scaled, slong prec)
{
    arb_t t, u, v;

    arb_init(t);
    arb_init(u);
    arb_init(v);

    arb_add_ui(t, nu, 1, prec);
    arb_hypgeom_bessel_k_integration(u, t, z, scaled, prec);
    arb_mul(u, u, t, prec);
    arb_mul_2exp_si(u, u, 1);
    arb_div(u, u, z, prec);

    arb_add_ui(t, nu, 2, prec);
    arb_hypgeom_bessel_k_integration(v, t, z, scaled, prec);
    arb_sub(res, v, u, prec);

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
}

void
arb_hypgeom_bessel_k_integration(arb_t res, const arb_t nu, const arb_t z, int scaled, slong prec)
{
    arb_t t, a, b, w;

    /* The integral for U(a, b, 2z) below only gives a finite enclosure
       for a >= 1 or 1 + a - b >= 1, that is, for |nu| >= 1/2. */
    if (arf_cmpabs_2exp_si(arb_midref(nu), -1) < 0)
    {
        bessel_k_integration_recurrence(res, nu, z, scaled, prec);
        return;
    }

    arb_init(t);
    arb_init(a);
    arb_init(b);
    arb_init(w);

    arb_one(a);
    arb_mul_2exp_si(a, a, -1);
    arb_add(a, a, nu, prec);

    arb_mul_2exp_si(b, nu, 1);
    arb_add_ui(b, b, 1, prec);

    arb_mul_2exp_si(w, z, 1);

    arb_hypgeom_u_integration(t, a, b, w, prec);

    if (arb_is_finite(t))
    {
        if (!scaled)
        {
            arb_neg(a, z);
            arb_exp(a, a, prec);
            arb_mul(t, t, a, prec);
        }

        arb_mul_2exp_si(w, z, 1);
        arb_pow(w, w, nu, prec);
        arb_mul(res, t, w, prec);

        arb_const_sqrt_pi(w, prec);
        arb_mul(res, res, w, prec);
    }
    else
    {
        arb_indeterminate(res);
    }

    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    arb_clear(w);
}
