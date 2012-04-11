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
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "padic_poly.h"

void padic_poly_realloc(padic_poly_t poly, long alloc, const fmpz_t p)
{
    if (alloc == 0)             /* Clear up, reinitialise */
    {
        padic_poly_clear(poly);
        padic_poly_init(poly);
        return;
    }

    if (poly->alloc)            /* Realloc */
    {
        padic_poly_truncate(poly, alloc, p);

        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, alloc * sizeof(fmpz));
        if (alloc > poly->alloc)
            mpn_zero((mp_ptr) (poly->coeffs + poly->alloc),
                     alloc - poly->alloc);
    }
    else                        /* Nothing allocated already so do it now */
    {
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    }

    poly->alloc = alloc;
}
