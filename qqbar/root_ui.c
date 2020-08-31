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
#include "qqbar.h"

void
qqbar_root_ui(qqbar_t res, const qqbar_t x, ulong n)
{
    if (n == 0)
    {
        flint_printf("qqbar_root_ui: n >= 1 is required");
        return;
    }
    else if (n == 1 || qqbar_is_zero(x) || qqbar_is_one(x))
    {
        qqbar_set(res, x);
    }
    else
    {
        slong i, d, prec, found;
        fmpz_poly_t H;
        fmpz_poly_factor_t fac;
        acb_t z, w, t;
        int pure_real;

        d = qqbar_degree(x);

        if (FLINT_BIT_COUNT(n) + FLINT_BIT_COUNT(d) > 30)
        {
            flint_printf("qqbar_root_ui: ludicrously high degree %wd * %wu", d, n);
            return;
        }

        /* handle principal roots of positive rational numbers */
        /* todo: could also handle conjugates of such roots */
        /* todo: also specialize roots of unity, etc. */
        if (_fmpz_vec_is_zero(QQBAR_COEFFS(x) + 1, d - 1) &&
            fmpz_sgn(QQBAR_COEFFS(x)) < 0 &&
            arb_contains_zero(acb_imagref(QQBAR_ENCLOSURE(x))) && arb_is_positive(acb_realref(QQBAR_ENCLOSURE(x))))
        {
            fmpq_t t;
            fmpq_init(t);
            fmpz_neg(fmpq_numref(t), QQBAR_COEFFS(x));
            fmpz_set(fmpq_denref(t), QQBAR_COEFFS(x) + d);
            qqbar_fmpq_root_ui(res, t, d * n);
            fmpq_clear(t);
            return;
        }

        fmpz_poly_init(H);
        fmpz_poly_factor_init(fac);
        acb_init(z);
        acb_init(w);
        acb_init(t);

        for (i = d; i >= 0; i--)
        {
            fmpz_poly_set_coeff_fmpz(H, i * n, QQBAR_COEFFS(x) + i);
        }

        fmpz_poly_factor(fac, H);
        acb_set(z, QQBAR_ENCLOSURE(x));
        pure_real = qqbar_is_real(x);

        for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, prec);
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
                if (_qqbar_validate_uniqueness(t, fac->p + found, w, 2 * prec))
                {
                    fmpz_poly_set(QQBAR_POLY(res), fac->p + found);
                    acb_set(QQBAR_ENCLOSURE(res), t);
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

