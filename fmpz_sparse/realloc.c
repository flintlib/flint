/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2008, 2009 William Hart
   
******************************************************************************/

#include "fmpz_sparse.h"

void
fmpz_sparse_realloc(fmpz_sparse_t poly, slong alloc)
{
    if (alloc == 0)             /* Clear up, reinitialise */
    {
    /*   fmpz_sparse_clear(poly);
        fmpz_sparse_init(poly);*/
        return;
    }

    if (poly->alloc)            /* Realloc */
    {
 /*       fmpz_sparse_truncate(poly, alloc);

        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, alloc * sizeof(fmpz));
        poly->expons = (fmpz *) flint_realloc(poly->coeffs, alloc * sizeof(fmpz));

        if (alloc > poly->alloc)
            flint_mpn_zero((mp_ptr) (poly->coeffs + poly->alloc),
                alloc - poly->alloc);
        if(alloc > poly->expons)
            flint_mpn_zero((mp_ptr) (poly->expons + poly->alloc),
                alloc - poly->alloc);*/
    }
    else                        /* Nothing allocated already so do it now */
    {
    /**    poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        poly->expons = (fmpz *) flint_calloc(alloc, sizeof(fmpz));*/
    }

    poly->alloc = alloc;
}
