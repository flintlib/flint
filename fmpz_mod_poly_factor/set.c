/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_factor_set(fmpz_mod_poly_factor_t res,
                    const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)
{
    if (res != fac)
    {
        if (fac->num == 0)
        {
            fmpz_mod_poly_factor_clear(res, ctx);
            fmpz_mod_poly_factor_init(res, ctx);
        }
        else
        {
            slong i;

            fmpz_mod_poly_factor_fit_length(res, fac->num, ctx);
            for (i = 0; i < fac->num; i++)
            {
                fmpz_mod_poly_set(res->poly + i, fac->poly + i, ctx);
                res->exp[i] = fac->exp[i];
            }
            for (; i < res->num; i++)
            {
                fmpz_mod_poly_zero(res->poly + i, ctx);
                res->exp[i] = 0;
            }
            res->num = fac->num;
        }
    }
}
