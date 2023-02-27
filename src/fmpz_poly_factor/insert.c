/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_poly.h"

void
fmpz_poly_factor_insert(fmpz_poly_factor_t fac, const fmpz_poly_t p, slong exp)
{
    slong i;

    for (i = 0; i < fac->num; i++)
    {
        if (fmpz_poly_equal(p, fac->p + i))
        {
            fac->exp[i] += exp;
            return;
        }
    }

    fmpz_poly_factor_fit_length(fac, i + 1);

    fmpz_poly_set(fac->p + i, p);
    fac->exp[i] = exp;
    fac->num = i + 1;
}
