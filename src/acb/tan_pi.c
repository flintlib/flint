/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_tan_pi(acb_t r, const acb_t z, slong prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_tan_pi(acb_realref(r), acb_realref(z), prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, prec + 4);
        arb_mul(t, acb_imagref(z), t, prec + 4);
        arb_tanh(acb_imagref(r), t, prec);
        arb_zero(acb_realref(r));
        arb_clear(t);
    }
    else
    {
        acb_t t;
        acb_init(t);

        if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) < 0)
        {
            acb_sin_cos_pi(r, t, z, prec + 4);
            acb_div(r, r, t, prec);
        }
        else
        {
            acb_mul_2exp_si(t, z, 1);

            if (arf_sgn(arb_midref(acb_imagref(z))) > 0)
            {
                acb_exp_pi_i(t, t, prec + 4);
                acb_add_ui(r, t, 1, prec + 4);
                acb_div(r, t, r, prec + 4);
                acb_mul_2exp_si(r, r, 1);
                acb_sub_ui(r, r, 1, prec);
                acb_div_onei(r, r);
            }
            else
            {
                acb_neg(t, t);
                acb_exp_pi_i(t, t, prec + 4);
                acb_add_ui(r, t, 1, prec + 4);
                acb_div(r, t, r, prec + 4);
                acb_mul_2exp_si(r, r, 1);
                acb_sub_ui(r, r, 1, prec);
                acb_mul_onei(r, r);
            }
        }

        acb_clear(t);
    }
}

