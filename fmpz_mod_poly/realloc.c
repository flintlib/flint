/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

void fmpz_mod_poly_realloc(fmpz_mod_poly_t poly, slong alloc,
                                                      const fmpz_mod_ctx_t ctx)
{
    if (alloc == 0)             /* Clear up, reinitialise */
    {
        if (poly->coeffs)
            _fmpz_vec_clear(poly->coeffs, poly->alloc);

        poly->coeffs = NULL;
        poly->length = 0;
        poly->alloc  = 0;

        return;
    }

    if (poly->alloc)            /* Realloc */
    {
        fmpz_mod_poly_truncate(poly, alloc, ctx);

        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, alloc * sizeof(fmpz));
        if (alloc > poly->alloc)
            flint_mpn_zero((mp_ptr) (poly->coeffs + poly->alloc),
                     alloc - poly->alloc);
    }
    else                        /* Nothing allocated already so do it now */
    {
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    }

    poly->alloc = alloc;
}
