/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"
#include "bernoulli.h"

static void
_arf_trunc(arf_t x)
{
    if (arf_sgn(x) < 0)
        arf_ceil(x, x);
    else
        arf_floor(x, x);
}

static void
acb_extract_bits(acb_t t, const acb_t z, slong b)
{
    acb_mul_2exp_si(t, z, b);
    _arf_trunc(arb_midref(acb_realref(t)));
    _arf_trunc(arb_midref(acb_imagref(t)));
    mag_zero(arb_radref(acb_realref(t)));
    mag_zero(arb_radref(acb_imagref(t)));
    acb_mul_2exp_si(t, t, -b);
}

void
acb_hypgeom_dilog_bitburst(acb_t res, acb_t z0, const acb_t z, slong prec)
{
    acb_t s, t, tprev, u;
    slong w;
    slong start = 16;

    acb_init(s);
    acb_init(t);
    acb_init(tprev);
    acb_init(u);

    acb_sub_ui(t, z, 1, 30);
    arb_abs(acb_imagref(t), acb_imagref(t));

    /* we don't want to end up on the branch cut */
    if (arb_contains_nonnegative(acb_realref(t))
        && !arb_gt(acb_imagref(t), acb_realref(t)))
    {
        acb_set(z0, z);
        acb_zero(res);
    }
    else
    {
        acb_extract_bits(t, z, start);
        acb_set(z0, t);
        acb_set(tprev, t);

        for (w = 2 * start; w < prec; w *= 2)
        {
            acb_extract_bits(t, z, w);
            acb_sub(u, t, z, 30);

            if (arf_cmpabs_2exp_si(arb_midref(acb_realref(u)), -prec / 8) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(u)), -prec / 8) < 0)
                break;

            acb_hypgeom_dilog_continuation(u, tprev, t, prec);
            acb_add(s, s, u, prec);
            acb_set(tprev, t);
        }

        acb_hypgeom_dilog_continuation(u, tprev, z, prec);
        acb_add(s, s, u, prec);

        acb_set(res, s);
    }

    acb_clear(s);
    acb_clear(t);
    acb_clear(tprev);
    acb_clear(u);
}

