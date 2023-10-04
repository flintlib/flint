/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_agm_rel_mag_err(arf_t m, arf_t eps, acb_srcptr a, arb_srcptr d,
    slong nb, slong prec)
{
    acb_t x, err;
    arb_t y;
    arf_t abs;
    slong k;

    acb_init(x);
    acb_init(err);
    arb_init(y);
    arf_init(abs);

    arf_zero(m);
    arf_zero(eps);

    for (k = 0; k < nb; k++)
    {
        arb_zero(y);
        arb_get_ubound_arf(arb_midref(y), &d[k], prec);
        arb_exp(y, y, prec);
        acb_mul_arb(x, &a[k], y, prec);

        acb_abs(y, x, prec);
        arb_get_ubound_arf(abs, y, prec);
        arf_max(m, m, abs);

        acb_zero(err);
        arf_set_mag(arb_midref(acb_realref(err)), arb_radref(acb_realref(x)));
        arf_set_mag(arb_midref(acb_imagref(err)), arb_radref(acb_imagref(x)));
        acb_abs(y, err, prec);
        arb_get_ubound_arf(abs, y, prec);
        arf_max(eps, eps, abs);
    }

    acb_clear(x);
    acb_clear(err);
    arb_clear(y);
    arf_clear(abs);
}
