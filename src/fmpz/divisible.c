/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2022, 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"

int
fmpz_divisible_ui(const fmpz_t f, ulong g)
{
    if (*f == WORD(0))
        return 1;

    if (!COEFF_IS_MPZ(*f))
    {
        ulong fa = FLINT_ABS(*f);
        return !(fa % g);
    }
    else
        return flint_mpz_divisible_ui_p(COEFF_TO_PTR(*f), g);
}

int fmpz_divisible(const fmpz_t x, const fmpz_t p)
{
    fmpz y = *x;
    fmpz q = *p;

    if (y == WORD(0))
        return 1;

    if (q == WORD(0))
        return 0;

    if (!COEFF_IS_MPZ(y))
    {
        if (!COEFF_IS_MPZ(q))
        {
            return !(y % q);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(q))
        {
            return flint_mpz_divisible_ui_p(COEFF_TO_PTR(y), FLINT_ABS(q));
        }
        else
        {
            return mpz_divisible_p(COEFF_TO_PTR(y), COEFF_TO_PTR(q));
        }
    }
}
