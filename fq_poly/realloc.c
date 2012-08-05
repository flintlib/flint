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

void
fq_poly_realloc(fq_poly_t poly, fq_ctx_t ctx, long alloc)
{

    if (alloc == 0)             /* Clear up, reinitialise */
    {
        fq_poly_clear(poly);
        fq_poly_init(poly,ctx);
        return;
    }

    if (poly->alloc)            /* Realloc */
    {
        fq_poly_truncate(poly, alloc);
        if( !fq_ctx_equal(poly->ctx,ctx)) fq_poly_change_ctx(poly,ctx);
        poly->coeffs = (fq_struct *) flint_realloc(poly->coeffs, alloc * sizeof(fq_t));
        if (alloc > poly->alloc)
            mpn_zero((mp_ptr) (poly->coeffs + poly->alloc),
                     alloc - poly->alloc);
    }
    else                        /* Nothing allocated already so do it now */
    {
        long i;
        if( !fq_ctx_equal(poly->ctx,ctx)) fq_poly_change_ctx(poly,ctx);
        poly->coeffs = (fq_struct *) flint_realloc(poly->coeffs, alloc * sizeof(fq_t));
        for(i=0;i<alloc;i++)
            fq_init2((poly->coeffs) + i,ctx);
    }

    poly->alloc = alloc;
}
