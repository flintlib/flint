/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

int
qqbar_pow(qqbar_t res, const qqbar_t x, const qqbar_t y)
{
    if (qqbar_is_zero(y))
    {
        qqbar_one(res);
        return 1;
    }
    else if (qqbar_is_one(y))
    {
        qqbar_set(res, x);
        return 1;
    }
    else if (qqbar_is_one(x))
    {
        qqbar_one(res);
        return 1;
    }
    else if (qqbar_is_zero(x))
    {
        int sign;

        sign = qqbar_sgn_re(y);

        if (sign > 0)
        {
            qqbar_zero(res);
            return 1;
        }

        return 0;
    }
    else if (qqbar_is_rational(y))
    {
        fmpq_t t;
        fmpz_t r;
        slong p;
        ulong q;
        int success;

        fmpq_init(t);
        fmpz_init(r);
        fmpz_neg(fmpq_numref(t), QQBAR_COEFFS(y));
        fmpz_set(fmpq_denref(t), QQBAR_COEFFS(y) + 1);

        success = 0;

        /* Fast path for roots of unity. */
        if (qqbar_is_root_of_unity(&p, &q, x))
        {
            fmpz_mul_si(fmpq_numref(t), fmpq_numref(t), p);
            fmpz_mul_ui(fmpq_denref(t), fmpq_denref(t), q);
            fmpz_mul_ui(r, fmpq_denref(t), 2);
            fmpz_fdiv_r(fmpq_numref(t), fmpq_numref(t), r);
            fmpq_canonicalise(t);

            if (COEFF_IS_MPZ(*fmpq_denref(t)))
            {
                flint_printf("qqbar_pow: excessive exponent\n");
                flint_abort();
            }

            qqbar_root_of_unity(res, *fmpq_numref(t), *fmpq_denref(t));
            success = 1;
        }

        if (!success)
        {
            if (COEFF_IS_MPZ(*fmpq_numref(t)) || COEFF_IS_MPZ(*fmpq_denref(t)))
            {
                flint_printf("qqbar_pow: excessive exponent\n");
                flint_abort();
            }

            p = *fmpq_numref(t);
            q = *fmpq_denref(t);

            qqbar_root_ui(res, x, q);
            if (p >= 0)
                qqbar_pow_ui(res, res, p);
            else
            {
                qqbar_pow_ui(res, res, -p);
                qqbar_inv(res, res);
            }

            success = 1;
        }

        fmpq_clear(t);
        fmpz_clear(r);

        return success;
    }
    else
    {
        return 0;
    }
}
