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

static int
check_root(fmpz_t p, const fmpz_t x, fmpz_t t, ulong d)
{
    if (fmpz_is_one(x))
    {
        fmpz_one(p);
        return 1;
    }
    else if (d == 2)
    {
        fmpz_sqrtrem(p, t, x);
        return fmpz_is_zero(t);
    }
    else
    {
        /* todo: need a rootrem function */
        fmpz_root(p, x, d);
        fmpz_pow_ui(t, p, d);
        return fmpz_equal(t, x);
    }
}

void
qqbar_fmpq_root_ui(qqbar_t res, const fmpq_t x, ulong b)
{
    ulong d;
    fmpz_t p, q, t;

    if (b == 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (b == 1 || fmpq_is_zero(x) || fmpq_is_one(x))
    {
        qqbar_set_fmpq(res, x);
        return;
    }

    if (b == 2)
    {
        if (fmpz_is_square(fmpq_numref(x)) && fmpz_is_square(fmpq_denref(x)))
        {
            fmpz_poly_fit_length(QQBAR_POLY(res), 2);
            _fmpz_poly_set_length(QQBAR_POLY(res), 2);
            fmpz_sqrt(QQBAR_COEFFS(res) + 0, fmpq_numref(x));
            fmpz_neg(QQBAR_COEFFS(res) + 0, QQBAR_COEFFS(res) + 0);
            fmpz_sqrt(QQBAR_COEFFS(res) + 1, fmpq_denref(x));
            acb_set_fmpz(QQBAR_ENCLOSURE(res), QQBAR_COEFFS(res) + 0);
            acb_neg(QQBAR_ENCLOSURE(res), QQBAR_ENCLOSURE(res));
            acb_div_fmpz(QQBAR_ENCLOSURE(res), QQBAR_ENCLOSURE(res), QQBAR_COEFFS(res) + 1, QQBAR_DEFAULT_PREC);
        }
        else
        {
            fmpz_poly_fit_length(QQBAR_POLY(res), 3);
            _fmpz_poly_set_length(QQBAR_POLY(res), 3);

            fmpz_set(QQBAR_COEFFS(res) + 0, fmpq_numref(x));
            fmpz_neg(QQBAR_COEFFS(res) + 0, QQBAR_COEFFS(res) + 0);
            fmpz_zero(QQBAR_COEFFS(res) + 1);
            fmpz_set(QQBAR_COEFFS(res) + 2, fmpq_denref(x));

            acb_set_fmpq(QQBAR_ENCLOSURE(res), x, QQBAR_DEFAULT_PREC);
            acb_sqrt(QQBAR_ENCLOSURE(res), QQBAR_ENCLOSURE(res), QQBAR_DEFAULT_PREC);
        }

        return;
    }

    if (fmpq_sgn(x) < 0)
    {
        qqbar_set_fmpq(res, x);
        qqbar_root_ui(res, res, b);
        return;
    }

    fmpz_init(p);
    fmpz_init(q);
    fmpz_init(t);

    for (d = b; d >= 1; d--)
    {
        if (d == 1)
        {
            fmpz_set(p, fmpq_numref(x));
            fmpz_set(q, fmpq_denref(x));
        }
        else if (b % d == 0)
        {
            if (check_root(p, fmpq_numref(x), t, d) && check_root(q, fmpq_denref(x), t, d))
            {
                b /= d;
                break;
            }
        }
    }

    fmpz_poly_set_fmpz(QQBAR_POLY(res), p);
    fmpz_poly_neg(QQBAR_POLY(res), QQBAR_POLY(res));
    fmpz_poly_set_coeff_fmpz(QQBAR_POLY(res), b, q);

    arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));
    arb_set_fmpz(acb_realref(QQBAR_ENCLOSURE(res)), p);
    arb_div_fmpz(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), q, QQBAR_DEFAULT_PREC);

    if (b != 1)
        arb_root_ui(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), b, QQBAR_DEFAULT_PREC);

    if (!arb_is_positive(acb_realref(QQBAR_ENCLOSURE(res))))
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(t);
}
