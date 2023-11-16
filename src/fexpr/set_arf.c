/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

void
fexpr_set_arf(fexpr_t res, const arf_t x)
{
    if (arf_is_zero(x))
    {
        fexpr_zero(res);
    }
    else if (arf_is_pos_inf(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_Infinity);
    }
    else if (arf_is_neg_inf(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_Infinity);
        fexpr_neg(res, res);
    }
    else if (arf_is_nan(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_Undefined);
    }
    else
    {
        fmpz_t m, e;
        fmpz_init(m);
        fmpz_init(e);
        arf_get_fmpz_2exp(m, e, x);

        if (0 <= *e && *e <= 20)
        {
            fmpz_mul_2exp(m, m, *e);
            fexpr_set_fmpz(res, m);
        }
        else if (-8 <= *e && *e < 0)
        {
            fmpq_t t;
            *fmpq_numref(t) = *m;
            *fmpq_denref(t) = (1 << (-*e));
            fexpr_set_fmpq(res, t);
        }
        else if (fmpz_is_pm1(m))
        {
            fexpr_t base, exp;
            fexpr_init(base);
            fexpr_init(exp);
            fexpr_set_si(base, 2);
            fexpr_set_fmpz(exp, e);
            fexpr_pow(res, base, exp);
            if (!fmpz_is_one(m))
                fexpr_neg(res, res);
            fexpr_clear(base);
            fexpr_clear(exp);
        }
        else
        {
            fexpr_t mantissa, base, exp;
            fexpr_init(mantissa);
            fexpr_init(base);
            fexpr_init(exp);
            fexpr_set_si(base, 2);
            fexpr_set_fmpz(exp, e);
            fexpr_pow(res, base, exp);
            fexpr_set_fmpz(mantissa, m);
            fexpr_mul(base, mantissa, res);
            fexpr_swap(res, base);
            fexpr_clear(mantissa);
            fexpr_clear(base);
            fexpr_clear(exp);
        }

        fmpz_clear(m);
        fmpz_clear(e);
    }
}
