/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz_poly.h"

void fmpz_poly_factor_set(fmpz_poly_factor_t res, const fmpz_poly_factor_t fac)
{
    if (res != fac)
    {
        if (fac->num == 0)
        {
            fmpz_poly_factor_clear(res);
            fmpz_poly_factor_init(res);
        }
        else
        {
            slong i;

            fmpz_poly_factor_fit_length(res, fac->num);
            fmpz_set(&(res->c), &(fac->c));
            for (i = 0; i < fac->num; i++)
            {
                fmpz_poly_set(res->p + i, fac->p + i);
                res->exp[i] = fac->exp[i];
            }
            for ( ; i < res->num; i++)
            {
                fmpz_poly_zero(res->p + i);
                res->exp[i] = 0;
            }
            res->num = fac->num;
        }
    }
}

