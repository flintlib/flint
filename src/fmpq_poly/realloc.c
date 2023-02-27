/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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
#include "fmpq_poly.h"

void fmpq_poly_realloc(fmpq_poly_t poly, slong alloc)
{
    if (alloc == 0)  /* Clear up, reinitialise */
    {
        fmpq_poly_clear(poly);
        fmpq_poly_init(poly);
        return;
    }

    if (poly->alloc)  /* Realloc */
    {
        if (poly->length > alloc)  /* Reduce the size */
        {
            slong i;
            for (i = alloc; i < poly->length; i++)
                _fmpz_demote(poly->coeffs + i);
            poly->length = alloc;
            _fmpq_poly_normalise(poly);
        }
        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, alloc * sizeof(fmpz));
        if (poly->alloc < alloc)
        {
            flint_mpn_zero((mp_ptr) (poly->coeffs + poly->alloc), alloc - poly->alloc);
        }
    }
    else  /* Nothing allocated, do it now */
    {
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    }
    
    poly->alloc = alloc;
}

