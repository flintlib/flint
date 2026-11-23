/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "arb_hypgeom/impl.h"
#include "acb_hypgeom/impl.h"

/* the working precision for the parameters and the gamma functions:
   with the bits of the largest parameter */
static slong
integration_wp(const arb_t a, const arb_t b, const arb_t z, slong prec)
{
    mag_t m, t;
    slong e;

    mag_init(m);
    mag_init(t);
    arb_get_mag(m, a);
    arb_get_mag(t, b);
    mag_max(m, m, t);
    arb_get_mag(t, z);
    mag_max(m, m, t);
    e = (mag_is_zero(m) || mag_is_inf(m)) ? 0 : MAG_EXP(m);
    mag_clear(m);
    mag_clear(t);
    return prec + 10 + FLINT_MIN(FLINT_MAX(e, 0), 4000);
}

/* The integral int_0^inf e^(-z t) t^a1 (1 + t)^ba1 dt: by the quadrature
   on the real interval (exp_log_integration.c), and for parameters
   which are too large for that, by the one for complex parameters. */
static int
_arb_hypgeom_u_integral_any(arb_t res, const arb_t a1, const arb_t ba1, const arb_t z, slong prec)
{
    acb_t I, ca, cb, cz;
    int ok;

    if (_arb_hypgeom_u_integral(res, a1, ba1, z, prec))
        return 1;

    acb_init(I); acb_init(ca); acb_init(cb); acb_init(cz);
    acb_set_arb(ca, a1);
    acb_set_arb(cb, ba1);
    acb_set_arb(cz, z);
    ok = _acb_hypgeom_u_integral(I, ca, cb, cz, prec);
    if (ok)
        arb_set(res, acb_realref(I));
    acb_clear(I); acb_clear(ca); acb_clear(cb); acb_clear(cz);
    return ok;
}

static int
_arb_hypgeom_u_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec)
{
    arb_t t, a1, ba1, I;
    slong wp = integration_wp(a, b, z, prec);
    int ok;

    arb_init(t);
    arb_init(a1);
    arb_init(ba1);
    arb_init(I);

    arb_sub_ui(a1, a, 1, wp);
    arb_sub(ba1, b, a, wp);
    arb_sub_ui(ba1, ba1, 1, wp);

    ok = arb_is_finite(z) && arb_is_positive(z);
    ok = ok && arb_is_nonnegative(a1);
    ok = ok && arb_is_finite(b);
    ok = ok && _arb_hypgeom_u_integral_any(I, a1, ba1, z, prec);

    if (!ok)
    {
        arb_indeterminate(res);
    }
    else
    {
        arb_rgamma(t, a, wp);
        arb_mul(res, I, t, wp);
        arb_set_round(res, res, prec);
    }

    arb_clear(t);
    arb_clear(a1);
    arb_clear(ba1);
    arb_clear(I);

    return ok;
}

void
arb_hypgeom_u_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec)
{
    arb_t res2;
    arb_init(res2);

    if (!_arb_hypgeom_u_integration(res2, a, b, z, prec))
    {
        arb_t c, d;
        arb_init(c);
        arb_init(d);
        arb_sub(c, a, b, prec);
        arb_add_ui(c, c, 1, prec);
        arb_sub_ui(d, b, 2, prec);
        arb_neg(d, d);

        if (_arb_hypgeom_u_integration(res2, c, d, z, prec))
        {
            arb_sub_ui(c, b, 1, prec);
            arb_neg(c, c);
            arb_pow(c, z, c, prec);
            arb_mul(res2, res2, c, prec);
        }

        arb_clear(c);
        arb_clear(d);
    }

    arb_swap(res, res2);
    arb_clear(res2);
}
