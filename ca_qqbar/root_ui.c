/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "ca_qqbar.h"

void
ca_qqbar_root_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong n)
{
    if (n == 0)
    {
        flint_printf("ca_qqbar_root_ui: n >= 1 is required");
        return;
    }
    else if (n == 1 || ca_qqbar_is_zero(x) || ca_qqbar_is_one(x))
    {
        ca_qqbar_set(res, x);
    }
    else
    {
        slong i, d, prec, found;
        fmpz_poly_t H;
        fmpz_poly_factor_t fac;
        acb_t z, w, t;
        int pure_real;

        d = ca_qqbar_degree(x);

        /* todo: fast handling of roots of rational numbers */

        if (FLINT_BIT_COUNT(n) + FLINT_BIT_COUNT(d) > 30)
        {
            flint_printf("ca_qqbar_root_ui: ludicrously high degree %wd * %wu", d, n);
            return;
        }

        fmpz_poly_init(H);
        fmpz_poly_factor_init(fac);
        acb_init(z);
        acb_init(w);
        acb_init(t);

        for (i = d; i >= 0; i--)
        {
            fmpz_poly_set_coeff_fmpz(H, i * n, CA_QQBAR_COEFFS(x) + i);
        }

        fmpz_poly_factor(fac, H);
        acb_set(z, CA_QQBAR_ENCLOSURE(x));
        pure_real = ca_qqbar_is_real(x);

        for (prec = CA_QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _ca_qqbar_enclosure_raw(z, CA_QQBAR_POLY(x), z, prec);
            if (pure_real)
                arb_zero(acb_imagref(z));

            acb_root_ui(w, z, n, prec);

            /* Look for potential roots -- we want exactly one */
            found = -1;
            for (i = 0; i < fac->num && found != -2; i++)
            {
                arb_fmpz_poly_evaluate_acb(t, fac->p + i, w, prec);
                if (acb_contains_zero(t))
                {
                    if (found == -1)
                        found = i;
                    else
                        found = -2;
                }
            }

            /* Check if the enclosure is good enough */
            if (found >= 0)
            {
                if (_ca_qqbar_validate_enclosure(t, fac->p + found, w, 2 * prec))
                {
                    fmpz_poly_set(CA_QQBAR_POLY(res), fac->p + found);
                    acb_set(CA_QQBAR_ENCLOSURE(res), t);
                    break;
                }
            }
        }

        fmpz_poly_clear(H);
        fmpz_poly_factor_clear(fac);
        acb_clear(z);
        acb_clear(w);
        acb_clear(t);
    }
}

