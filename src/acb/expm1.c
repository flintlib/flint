/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_expm1(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_real(z))
    {
        arb_expm1(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
    }
    else if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -3) <= 0 &&
             arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -3) <= 0 &&
             arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -3) <= 0 &&
             arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -3) <= 0)
    {
        arf_srcptr midmax;
        slong extra;

        if (arf_cmpabs(arb_midref(acb_realref(z)), arb_midref(acb_imagref(z))) >= 0)
            midmax = arb_midref(acb_realref(z));
        else
            midmax = arb_midref(acb_imagref(z));

        if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -prec - 100) > 0)
        {
            extra = -ARF_EXP(midmax);
            extra = FLINT_MIN(extra, prec + 100);
            extra = FLINT_MAX(extra, 0);
            acb_exp(res, z, prec + extra + 4);
            acb_sub_ui(res, res, 1, prec);
        }
        else
        {
            /* lazy solution: e^z-1 = 4 (sinh(z/4)+cosh(z/4))^2 sinh(z/4) cosh(z/4) */
            acb_t t, u;
            acb_init(t);
            acb_init(u);
            acb_mul_2exp_si(t, z, -2);
            acb_sinh_cosh(t, u, t, prec + 4);
            acb_add(res, t, u, prec + 4);
            acb_mul(res, res, res, prec + 4);
            acb_mul(t, t, u, prec + 4);
            acb_mul(res, res, t, prec);
            acb_mul_2exp_si(res, res, 2);
            acb_clear(t);
            acb_clear(u);
        }
    }
    else
    {
        acb_exp(res, z, prec + 4);
        acb_sub_ui(res, res, 1, prec);
    }
}

