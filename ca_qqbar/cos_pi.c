/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_cos_pi(ca_qqbar_t res, slong p, ulong q)
{
    fmpq_t t;
    ulong a, b;
    slong prec;

    fmpq_init(t);

    if (q == 0)
    {
        flint_printf("ca_qqbar_cos_pi: q = 0\n");
        flint_abort();
    }

    fmpq_set_si(t, p, q);
    fmpq_div_2exp(t, t, 1);
    fmpz_fdiv_r(fmpq_numref(t), fmpq_numref(t), fmpq_denref(t));

    a = fmpz_get_ui(fmpq_numref(t));
    b = fmpz_get_ui(fmpq_denref(t));

    if (a == 0)
    {
        ca_qqbar_one(res);
    }
    else if (b == 2)
    {
        ca_qqbar_set_si(res, -1);
    }
    else if (b == 3)
    {
        ca_qqbar_one(res);
        ca_qqbar_neg(res, res);
        ca_qqbar_mul_2exp_si(res, res, -1);
    }
    else if (b == 4)
    {
        ca_qqbar_zero(res);
    }
    else if (b == 6)
    {
        ca_qqbar_one(res);
        ca_qqbar_mul_2exp_si(res, res, -1);
    }
    else
    {
        fmpz_poly_cos_minpoly(CA_QQBAR_POLY(res), b);
        fmpq_mul_2exp(t, t, 1);

        for (prec = CA_QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            arb_cos_pi_fmpq(acb_realref(CA_QQBAR_ENCLOSURE(res)), t, prec);
            arb_zero(acb_imagref(CA_QQBAR_ENCLOSURE(res)));
            acb_mul_2exp_si(CA_QQBAR_ENCLOSURE(res), CA_QQBAR_ENCLOSURE(res), 1);

            if (_ca_qqbar_validate_uniqueness(CA_QQBAR_ENCLOSURE(res),
                    CA_QQBAR_POLY(res), CA_QQBAR_ENCLOSURE(res), prec * 2))
            {
                break;
            }
        }

        ca_qqbar_mul_2exp_si(res, res, -1);
    }

    fmpq_clear(t);
}

