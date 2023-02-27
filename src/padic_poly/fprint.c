/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>

#include "padic_poly.h"

int _padic_poly_fprint(FILE *file, const fmpz *poly, slong val, slong len, 
                       const padic_ctx_t ctx)
{
    slong i, v;
    fmpz_t u;

    if (len == 0)
    {
        flint_fprintf(file, "0");
        return 1;
    }

    fmpz_init(u);

    flint_fprintf(file, "%wd ", len);

    for (i = 0; i < len; i++)
    {
        flint_fprintf(file, " ");

        if (fmpz_is_zero(poly + i))
        {
            flint_fprintf(file, "0");
        }
        else
        {
            v = val + fmpz_remove(u, poly + i, ctx->p);

            _padic_fprint(file, u, v, ctx);
        }
    }

    fmpz_clear(u);

    return 1;
}

int padic_poly_fprint(FILE *file, const padic_poly_t poly, 
                      const padic_ctx_t ctx)
{
    _padic_poly_fprint(file, poly->coeffs, poly->val, poly->length, ctx);

    return 1;
}

