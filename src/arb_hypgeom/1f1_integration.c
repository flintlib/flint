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

/* The integral int_0^1 e^(z t) t^a1 (1 - t)^ba1 dt: by the quadrature on
   the real interval (exp_log_integration.c), and for parameters which
   are too large for that, by the one for complex parameters
   (acb_hypgeom/exp_log_integration.c), which works relative to the
   saddle point. */
static int
_arb_hypgeom_1f1_integral_any(arb_t res, const arb_t a1, const arb_t ba1, const arb_t z, slong prec)
{
    acb_t I, ca, cb, cz;
    int ok;

    if (_arb_hypgeom_1f1_integral(res, a1, ba1, z, prec))
        return 1;

    acb_init(I); acb_init(ca); acb_init(cb); acb_init(cz);
    acb_set_arb(ca, a1);
    acb_set_arb(cb, ba1);
    acb_set_arb(cz, z);
    ok = _acb_hypgeom_1f1_integral(I, ca, cb, cz, prec);
    if (ok)
        arb_set(res, acb_realref(I));
    acb_clear(I); acb_clear(ca); acb_clear(cb); acb_clear(cz);
    return ok;
}

void
arb_hypgeom_1f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)
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

    ok = arb_is_finite(z);
    ok = ok && arb_is_nonnegative(a1);
    ok = ok && arb_is_nonnegative(ba1);
    ok = ok && _arb_hypgeom_1f1_integral_any(I, a1, ba1, z, prec);

    if (!ok)
    {
        arb_indeterminate(res);
    }
    else
    {
        if (!regularized)
        {
            arb_gamma(t, b, wp);
            arb_mul(I, I, t, wp);
        }
        arb_rgamma(t, a, wp);
        arb_mul(I, I, t, wp);
        arb_sub(t, b, a, wp);
        arb_rgamma(t, t, wp);
        arb_mul(res, I, t, wp);
        arb_set_round(res, res, prec);
    }

    arb_clear(t);
    arb_clear(a1);
    arb_clear(ba1);
    arb_clear(I);
}
