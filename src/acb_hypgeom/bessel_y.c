/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

static void
phase(acb_t res, const arb_t re, const arb_t im)
{
    if (arb_is_nonnegative(re) || arb_is_negative(im))
    {
        acb_one(res);
    }
    else if (arb_is_negative(re) && arb_is_nonnegative(im))
    {
        acb_set_si(res, -3);
    }
    else
    {
        arb_zero(acb_imagref(res));
        /* -1 +/- 2 */
        arf_set_si(arb_midref(acb_realref(res)), -1);
        mag_one(arb_radref(acb_realref(res)));
        mag_mul_2exp_si(arb_radref(acb_realref(res)), arb_radref(acb_realref(res)), 1);
    }
}

void
acb_hypgeom_bessel_jy(acb_t res1, acb_t res2, const acb_t nu, const acb_t z, slong prec)
{
    acb_t jnu, t, u, v;

    /* both j and y are zero at -inf and +inf when nu is finite */
    if (acb_is_finite(nu) && !acb_is_finite(z) &&
        acb_is_real(z) && !acb_contains_zero(z))
    {
        if (res1 != NULL) acb_zero(res1);
        acb_zero(res2);
        return;
    }

    acb_init(jnu);
    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_hypgeom_bessel_j(jnu, nu, z, prec);

    if (acb_is_int(nu))
    {
        int is_real = acb_is_real(nu) && acb_is_real(z)
            && arb_is_positive(acb_realref(z));

        acb_mul_onei(t, z);
        acb_hypgeom_bessel_k(t, nu, t, prec);
        acb_onei(u);
        acb_pow(u, u, nu, prec);
        acb_mul(t, t, u, prec);
        acb_const_pi(u, prec);
        acb_div(t, t, u, prec);
        acb_mul_2exp_si(t, t, 1);
        acb_neg(t, t);

        phase(v, acb_realref(z), acb_imagref(z));
        acb_mul(u, jnu, v, prec);
        acb_mul_onei(u, u);

        acb_sub(res2, t, u, prec);

        if (is_real)
            arb_zero(acb_imagref(res2));
    }
    else
    {
        acb_sin_cos_pi(t, u, nu, prec);
        acb_mul(v, jnu, u, prec);
        acb_neg(u, nu);
        acb_hypgeom_bessel_j(u, u, z, prec);
        acb_sub(v, v, u, prec);
        acb_div(res2, v, t, prec);
    }

    if (res1 != NULL)
        acb_set(res1, jnu);

    acb_clear(jnu);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_hypgeom_bessel_y(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    acb_hypgeom_bessel_jy(NULL, res, nu, z, prec);
}

