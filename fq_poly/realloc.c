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
    Copyright (C) 2012 Andres Goens
   
******************************************************************************/

#include "fq_poly.h"

void fq_poly_realloc(fq_poly_t poly, long alloc)
{
    long i;

    if (alloc == 0)             /* Clear up, reinitialise */
    {
        fq_poly_clear(poly);
        fq_poly_init(poly);
    }
    else if (poly->alloc)       /* Realloc */
    {
        for (i = alloc; i < poly->alloc; i++)
            fq_clear(poly->coeffs + i);

        poly->coeffs = (fq_struct *) flint_realloc(poly->coeffs, alloc * sizeof(fq_struct));

        for (i = poly->alloc; i < alloc; i++)
            fq_init(poly->coeffs + i);

        poly->length = FLINT_MIN(poly->length, alloc);
        _fq_poly_normalise(poly);
    }
    else                        /* Nothing allocated already so do it now */
    {
        poly->coeffs = (fq_struct *) flint_malloc(alloc * sizeof(fq_struct));

        for(i = 0; i < alloc; i++)
            fq_init(poly->coeffs + i);
    }
    poly->alloc = alloc;
}
