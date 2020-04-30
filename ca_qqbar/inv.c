/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "ca_qqbar.h"

void
ca_qqbar_inv(ca_qqbar_t res, const ca_qqbar_t x)
{
    slong d;

    if (ca_qqbar_is_zero(x))
    {
        flint_printf("ca_qqbar_inv: division by zero\n");
        flint_abort();
    }

    if (ca_qqbar_is_one(x) || ca_qqbar_is_neg_one(x))
    {
        ca_qqbar_set(res, x);
        return;
    }

    d = ca_qqbar_degree(x);

    if (d == 1)
    {
        fmpz_poly_reverse(CA_QQBAR_POLY(res), CA_QQBAR_POLY(x), d + 1);
        if (fmpz_sgn(CA_QQBAR_COEFFS(res) + d) < 0)
            fmpz_poly_neg(CA_QQBAR_POLY(res), CA_QQBAR_POLY(res));

        arb_fmpz_div_fmpz(acb_realref(CA_QQBAR_ENCLOSURE(res)),
            CA_QQBAR_COEFFS(res), CA_QQBAR_COEFFS(res) + 1, CA_QQBAR_DEFAULT_PREC);
        arb_neg(acb_realref(CA_QQBAR_ENCLOSURE(res)), acb_realref(CA_QQBAR_ENCLOSURE(res)));
        arb_zero(acb_imagref(CA_QQBAR_ENCLOSURE(res)));
    }
    else
    {
        fmpz_poly_t pol;
        slong prec;
        acb_t z, t;

        fmpz_poly_init(pol);
        acb_init(z);
        acb_init(t);

        fmpz_poly_reverse(pol, CA_QQBAR_POLY(x), d + 1);
        if (fmpz_sgn(pol->coeffs + d) < 0)
            fmpz_poly_neg(pol, pol);
        acb_set(z, CA_QQBAR_ENCLOSURE(x));

        for (prec = CA_QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _ca_qqbar_enclosure_raw(z, CA_QQBAR_POLY(x), z, prec);
            acb_inv(t, z, prec);

            if (_ca_qqbar_validate_enclosure(t, pol, t, 2 * prec))
            {
                fmpz_poly_set(CA_QQBAR_POLY(res), pol);
                acb_set(CA_QQBAR_ENCLOSURE(res), t);
                break;
            }

        }

        fmpz_poly_clear(pol);
        acb_clear(z);
        acb_clear(t);
    }
}

