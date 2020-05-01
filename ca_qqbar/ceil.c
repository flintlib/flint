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
ca_qqbar_ceil(fmpz_t res, const ca_qqbar_t x)
{
    if (ca_qqbar_is_rational(x))
    {
        fmpz_fdiv_q(res, CA_QQBAR_COEFFS(x), CA_QQBAR_COEFFS(x) + 1);
        fmpz_neg(res, res);
    }
    else
    {
        arb_t v;
        arb_init(v);

        arb_ceil(v, acb_realref(CA_QQBAR_ENCLOSURE(x)), CA_QQBAR_DEFAULT_PREC);
        if (!arb_get_unique_fmpz(res, v))
        {
            mag_t t;
            slong size, prec;
            acb_t z;

            mag_init(t);
            acb_init(z);

            acb_get_mag(t, CA_QQBAR_ENCLOSURE(x));
            if (mag_cmp_2exp_si(t, 0) < 0)
                mag_one(t);
            size = *MAG_EXPREF(t);
            prec = FLINT_MAX(CA_QQBAR_DEFAULT_PREC * 2, 2 * size + 32);

            acb_set(z, CA_QQBAR_ENCLOSURE(x));
            _ca_qqbar_enclosure_raw(z, CA_QQBAR_POLY(x), z, prec);
            arb_ceil(v, acb_realref(z), prec);

            /* Do an exact computation */
            if (!arb_get_unique_fmpz(res, v))
            {
                ca_qqbar_t u;
                ca_qqbar_init(u);

                arb_set_d(v, -0.5);
                arb_add(v, v, acb_realref(z), prec);
                arb_ceil(v, v, prec);

                if (!arb_get_unique_fmpz(res, v))
                {
                    flint_printf("ca_qqbar_ceil: either ceil(x) or ceil(x-1/2) should evaluate numerically\n");
                    flint_abort();
                }

                ca_qqbar_set_fmpz(u, res);
                ca_qqbar_sub(u, x, u);

                if (ca_qqbar_sgn_re(u) > 0)
                    fmpz_add_ui(res, res, 1);

                ca_qqbar_clear(u);
            }

            mag_clear(t);
            acb_clear(z);
        }

        arb_clear(v);
    }
}

