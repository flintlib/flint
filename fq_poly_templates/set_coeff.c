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
    Copyright (C) 2013 Mike Hansen
   
******************************************************************************/


#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_set_coeff) (TEMPLATE(T, poly_t) poly, slong n,
                             const TEMPLATE(T, t) x,
                             const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    TEMPLATE(T, poly_fit_length) (poly, n + 1, ctx);

    if (n + 1 > poly->length)   /* Insert zeros if needed */
    {
        for (i = poly->length; i < n; i++)
            TEMPLATE(T, zero) (poly->coeffs + i, ctx);
        poly->length = n + 1;
    }

    TEMPLATE(T, set) (poly->coeffs + n, x, ctx);
    _TEMPLATE(T, poly_normalise) (poly, ctx);
}


#endif
