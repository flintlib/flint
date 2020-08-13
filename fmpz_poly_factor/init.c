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

void fmpz_poly_factor_init(fmpz_poly_factor_t fac)
{
    fmpz_init_set_ui(&(fac->c), 1);
    fac->p     = NULL;
    fac->exp   = NULL;
    fac->num   = 0;
    fac->alloc = 0;
}

void fmpz_poly_factor_init2(fmpz_poly_factor_t fac, slong alloc)
{
    fmpz_init_set_ui(&(fac->c), 1);

    if (alloc)
    {
        slong i;

        fac->p   = flint_malloc(alloc * sizeof(fmpz_poly_struct));
        fac->exp = flint_malloc(alloc * sizeof(slong));

        for (i = 0; i < alloc; i++)
        {
            fmpz_poly_init(fac->p + i);
            fac->exp[i] = WORD(0);
        }
    }
    else
    {
        fac->p   = NULL;
        fac->exp = NULL;
    }

    fac->num   = 0;
    fac->alloc = alloc;
}

