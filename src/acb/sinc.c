/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

static int
_acb_is_large(const acb_t z)
{
    int res;
    mag_t t;
    mag_init(t);
    acb_get_mag_lower(t, z);
    res = mag_cmp_2exp_si(t, 0) >= 0;
    mag_clear(t);
    return res;
}

static void
_acb_sinc_direct(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_zero(z))
    {
        acb_one(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_sin(t, z, prec + 2);
        acb_div(res, t, z, prec);
        acb_clear(t);
    }
}

static void
_acb_sinc_diffbound(acb_t res, const acb_t z, slong prec)
{
    mag_t u, v;
    int isreal;

    mag_init(u);
    mag_init(v);

    isreal = arb_is_zero(acb_realref(z));

    /* exp(|im(z)|) is a crude bound for |sinc'(z)| */
    arb_get_mag(u, acb_imagref(z));
    mag_hypot(v, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
    mag_exp(u, u);
    mag_mul(u, u, v);

    arf_set(arb_midref(acb_realref(res)), arb_midref(acb_realref(z)));
    arf_set(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(z)));
    mag_zero(arb_radref(acb_realref(res)));
    mag_zero(arb_radref(acb_imagref(res)));

    _acb_sinc_direct(res, res, prec);

    if (isreal)
        arb_add_error_mag(acb_realref(res), u);
    else
        acb_add_error_mag(res, u);

    mag_clear(u);
    mag_clear(v);
}

void
acb_sinc(acb_t res, const acb_t z, slong prec)
{
    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
    }
    else if (acb_is_real(z))
    {
        arb_sinc(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
    }
    else if (acb_is_exact(z) || _acb_is_large(z))
    {
        _acb_sinc_direct(res, z, prec);
    }
    else
    {
        _acb_sinc_diffbound(res, z, prec);
    }
}

