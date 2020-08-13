/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
fmpz_cmp(const fmpz_t f, const fmpz_t g)
{
    int sign;

    if (f == g)
        return 0;  /* aliased inputs */

    if (!COEFF_IS_MPZ(*f))
    {
        if (!COEFF_IS_MPZ(*g))
        {
            return (*f < *g ? -1 : *f > *g);
        }
        else  /* f is small, g is large */
        {
            sign = mpz_sgn(COEFF_TO_PTR(*g));
            return (*f >= 0 && sign < 0) ? 1 : -sign;
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(*g))  /* f is large, and g is small */
        {
            sign = mpz_sgn(COEFF_TO_PTR(*f));
            return (*g >= 0 && sign < 0) ? -1 : sign;
        }
        else
            return mpz_cmp(COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
    }
}
