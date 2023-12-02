/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_floor(fmpz_t res, const qqbar_t x)
{
    if (qqbar_is_rational(x))
    {
        fmpz_cdiv_q(res, QQBAR_COEFFS(x), QQBAR_COEFFS(x) + 1);
        fmpz_neg(res, res);
    }
    else
    {
        arb_t v;
        arb_init(v);

        arb_floor(v, acb_realref(QQBAR_ENCLOSURE(x)), QQBAR_DEFAULT_PREC);
        if (!arb_get_unique_fmpz(res, v))
        {
            mag_t t;
            slong size, prec;
            acb_t z;

            mag_init(t);
            acb_init(z);

            acb_get_mag(t, QQBAR_ENCLOSURE(x));
            if (mag_cmp_2exp_si(t, 0) < 0)
                mag_one(t);
            size = *MAG_EXPREF(t);
            prec = FLINT_MAX(QQBAR_DEFAULT_PREC * 2, 2 * size + 32);

            acb_set(z, QQBAR_ENCLOSURE(x));
            _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, prec);
            arb_floor(v, acb_realref(z), prec);

            /* Do an exact computation */
            if (!arb_get_unique_fmpz(res, v))
            {
                qqbar_t u;
                qqbar_init(u);

                arb_set_d(v, 0.5);
                arb_add(v, v, acb_realref(z), prec);
                arb_floor(v, v, prec);

                if (!arb_get_unique_fmpz(res, v))
                {
                    flint_throw(FLINT_ERROR, "qqbar_floor: either floor(x) or floor(x+1/2) should evaluate numerically\n");
                }

                qqbar_set_fmpz(u, res);
                qqbar_sub(u, x, u);

                if (qqbar_sgn_re(u) < 0)
                    fmpz_sub_ui(res, res, 1);

                qqbar_clear(u);
            }

            mag_clear(t);
            acb_clear(z);
        }

        arb_clear(v);
    }
}

