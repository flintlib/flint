/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>

#include "fmpz.h"
#include "fmpz_mod_poly.h"

int _fmpz_mod_poly_fprint(FILE * file, const fmpz *poly, slong len, 
                          const fmpz_t p)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%wd ", len);
    if (r <= 0)
        return r;

    r = fmpz_fprint(file, p);
    if (r <= 0)
        return r;

    if (len == 0)
        return r;

    r = flint_fprintf(file, " ");
    if (r <= 0)
        return r;

    for (i = 0; (r > 0) && (i < len); i++)
    {
        r = flint_fprintf(file, " ");
        if (r <= 0)
            return r;
        r = fmpz_fprint(file, poly + i);
        if (r <= 0)
            return r;
    }

    return r;
}

int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly,
                                                      const fmpz_mod_ctx_t ctx)
{
   return _fmpz_mod_poly_fprint(file, poly->coeffs, poly->length,
                                                    fmpz_mod_ctx_modulus(ctx));
}

