/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_poly.h"
#include "qqbar.h"

void
qqbar_cos_pi(qqbar_t res, slong p, ulong q)
{
    fmpq_t t;
    ulong a, b;
    slong prec;

    fmpq_init(t);

    if (q == 0)
    {
        flint_throw(FLINT_ERROR, "qqbar_cos_pi: q = 0\n");
    }

    fmpq_set_si(t, p, q);
    fmpq_div_2exp(t, t, 1);
    fmpz_fdiv_r(fmpq_numref(t), fmpq_numref(t), fmpq_denref(t));

    a = fmpz_get_ui(fmpq_numref(t));
    b = fmpz_get_ui(fmpq_denref(t));

    if (a == 0)
    {
        qqbar_one(res);
    }
    else if (b == 2)
    {
        qqbar_set_si(res, -1);
    }
    else if (b == 3)
    {
        qqbar_one(res);
        qqbar_neg(res, res);
        qqbar_mul_2exp_si(res, res, -1);
    }
    else if (b == 4)
    {
        qqbar_zero(res);
    }
    else if (b == 6)
    {
        qqbar_one(res);
        qqbar_mul_2exp_si(res, res, -1);
    }
    else
    {
        fmpz_poly_cos_minpoly(QQBAR_POLY(res), b);
        fmpq_mul_2exp(t, t, 1);

        for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            arb_cos_pi_fmpq(acb_realref(QQBAR_ENCLOSURE(res)), t, prec);
            arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));
            acb_mul_2exp_si(QQBAR_ENCLOSURE(res), QQBAR_ENCLOSURE(res), 1);

            if (_qqbar_validate_uniqueness(QQBAR_ENCLOSURE(res),
                    QQBAR_POLY(res), QQBAR_ENCLOSURE(res), prec * 2))
            {
                break;
            }
        }

        qqbar_mul_2exp_si(res, res, -1);
    }

    fmpq_clear(t);
}

