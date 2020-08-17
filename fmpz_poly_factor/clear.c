/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_factor_clear(fmpz_poly_factor_t fac)
{
    if (fac->alloc)
    {
        slong i;

        for (i = 0; i < fac->alloc; i++)
        {
            fmpz_poly_clear(fac->p + i);
        }

        flint_free(fac->p);
        flint_free(fac->exp);
        fac->p   = NULL;
        fac->exp = NULL;
    }

    fmpz_clear(&(fac->c));
}

