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
integration_wp(const arb_t a, const arb_t b, const arb_t c, const arb_t z, slong prec)
{
    mag_t m, t;
    slong e;

    mag_init(m);
    mag_init(t);
    arb_get_mag(m, a);
    arb_get_mag(t, b);
    mag_max(m, m, t);
    arb_get_mag(t, c);
    mag_max(m, m, t);
    arb_get_mag(t, z);
    mag_max(m, m, t);
    e = (mag_is_zero(m) || mag_is_inf(m)) ? 0 : MAG_EXP(m);
    mag_clear(m);
    mag_clear(t);
    return prec + 10 + FLINT_MIN(FLINT_MAX(e, 0), 4000);
}

/* The integral int_0^1 t^b1 (1 - t)^cb1 (1 - z t)^nega dt: by the
   quadrature on the real interval (exp_log_integration.c), and for
   parameters which are too large for that, by the one for complex
   parameters. */
static int
_arb_hypgeom_2f1_integral_any(arb_t res, const arb_t b1, const arb_t cb1, const arb_t nega,
    const arb_t z, slong prec)
{
    acb_t I, cb, cc, ca, cz;
    int ok;

    if (_arb_hypgeom_2f1_integral(res, b1, cb1, nega, z, prec))
        return 1;

    acb_init(I); acb_init(cb); acb_init(cc); acb_init(ca); acb_init(cz);
    acb_set_arb(cb, b1);
    acb_set_arb(cc, cb1);
    acb_set_arb(ca, nega);
    acb_set_arb(cz, z);
    ok = _acb_hypgeom_2f1_integral(I, cb, cc, ca, cz, prec);
    if (ok)
        arb_set(res, acb_realref(I));
    acb_clear(I); acb_clear(cb); acb_clear(cc); acb_clear(ca); acb_clear(cz);
    return ok;
}

static int
_arb_hypgeom_2f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec)
{
    arb_t t, b1, cb1, nega, I;
    slong wp = integration_wp(a, b, c, z, prec);
    int ok;

    arb_init(t);
    arb_init(b1);
    arb_init(cb1);
    arb_init(nega);
    arb_init(I);

    arb_sub_ui(b1, b, 1, wp);
    arb_sub(cb1, c, b, wp);
    arb_sub_ui(cb1, cb1, 1, wp);
    arb_neg(nega, a);

    arb_one(t);
    ok = arb_is_finite(z) && arb_lt(z, t);
    ok = ok && arb_is_nonnegative(b1);
    ok = ok && arb_is_nonnegative(cb1);
    ok = ok && _arb_hypgeom_2f1_integral_any(I, b1, cb1, nega, z, prec);

    if (!ok)
    {
        arb_indeterminate(res);
    }
    else
    {
        if (!(regularized & 1))
        {
            arb_gamma(t, c, wp);
            arb_mul(I, I, t, wp);
        }
        arb_rgamma(t, b, wp);
        arb_mul(I, I, t, wp);
        arb_sub(t, c, b, wp);
        arb_rgamma(t, t, wp);
        arb_mul(res, I, t, wp);
        arb_set_round(res, res, prec);
    }

    arb_clear(t);
    arb_clear(b1);
    arb_clear(cb1);
    arb_clear(nega);
    arb_clear(I);

    return ok;
}

void
arb_hypgeom_2f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec)
{
    arb_t res2;
    arb_init(res2);

    if (arf_cmp(arb_midref(a), arb_midref(b)) < 0)
    {
        if (!_arb_hypgeom_2f1_integration(res2, a, b, c, z, regularized, prec))
            _arb_hypgeom_2f1_integration(res2, b, a, c, z, regularized, prec);
    }
    else
    {
        if (!_arb_hypgeom_2f1_integration(res2, b, a, c, z, regularized, prec))
            _arb_hypgeom_2f1_integration(res2, a, b, c, z, regularized, prec);
    }

    arb_swap(res, res2);
    arb_clear(res2);
}
