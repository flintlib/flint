/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

void
fmpz_factor_divisor_sigma(fmpz_t res, ulong k, const fmpz_factor_t fac)
{
    slong i;

    fmpz_one(res);

    if (fac->num == 0)
        return;

    if (k == 0)
    {
        for (i = 0; i < fac->num; i++)
            fmpz_mul_ui(res, res, fac->exp[i] + 1);
    }
    else
    {
        fmpz * p;
        fmpz_t r;

        p = _fmpz_vec_init(fac->num);
        fmpz_init(r);

        for (i = 0; i < fac->num; i++)
        {
            fmpz_pow_ui(p + i, fac->p + i, k);
            fmpz_pow_ui(r, p + i, fac->exp[i]  + 1);
            fmpz_sub_ui(r, r, 1);
            fmpz_sub_ui(p + i, p + i, 1);
            fmpz_divexact(p + i, r, p + i);
        }

        _fmpz_vec_prod(res, p, fac->num);

        _fmpz_vec_clear(p, fac->num);
        fmpz_clear(r);
    }
}

void
fmpz_divisor_sigma(fmpz_t res, ulong k, const fmpz_t n)
{
    fmpz_factor_t fac;

    if (fmpz_is_zero(n))
    {
        fmpz_zero(res);
        return;
    }

    fmpz_factor_init(fac);
    fmpz_factor(fac, n);
    fmpz_factor_divisor_sigma(res, k, fac);
    fmpz_factor_clear(fac);
}

