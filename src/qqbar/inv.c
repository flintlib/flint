/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "qqbar.h"

void
qqbar_inv(qqbar_t res, const qqbar_t x)
{
    slong d;

    if (qqbar_is_zero(x))
    {
        flint_throw(FLINT_ERROR, "qqbar_inv: division by zero\n");
    }

    if (qqbar_is_one(x) || qqbar_is_neg_one(x))
    {
        qqbar_set(res, x);
        return;
    }

    d = qqbar_degree(x);

    if (d == 1)
    {
        fmpz_poly_reverse(QQBAR_POLY(res), QQBAR_POLY(x), d + 1);
        if (fmpz_sgn(QQBAR_COEFFS(res) + d) < 0)
            fmpz_poly_neg(QQBAR_POLY(res), QQBAR_POLY(res));

        arb_fmpz_div_fmpz(acb_realref(QQBAR_ENCLOSURE(res)),
            QQBAR_COEFFS(res), QQBAR_COEFFS(res) + 1, QQBAR_DEFAULT_PREC);
        arb_neg(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)));
        arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));
    }
    else
    {
        fmpz_poly_t pol;
        slong prec;
        acb_t z, t;

        fmpz_poly_init(pol);
        acb_init(z);
        acb_init(t);

        fmpz_poly_reverse(pol, QQBAR_POLY(x), d + 1);
        if (fmpz_sgn(pol->coeffs + d) < 0)
            fmpz_poly_neg(pol, pol);
        acb_set(z, QQBAR_ENCLOSURE(x));

        for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, prec);
            acb_inv(t, z, prec);

            if (_qqbar_validate_uniqueness(t, pol, t, 2 * prec))
            {
                fmpz_poly_set(QQBAR_POLY(res), pol);
                acb_set(QQBAR_ENCLOSURE(res), t);
                break;
            }

        }

        fmpz_poly_clear(pol);
        acb_clear(z);
        acb_clear(t);
    }
}

