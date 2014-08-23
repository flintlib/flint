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

    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
   
******************************************************************************/


#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_set_trunc) (TEMPLATE(T, poly_t) poly1, TEMPLATE(T, poly_t) poly2, slong len,
                            const TEMPLATE(T, ctx_t) ctx)
{
    if (poly1 == poly2)
        TEMPLATE(T, poly_truncate) (poly1, len, ctx);
    else if (len >= poly2->length)
        TEMPLATE(T, poly_set) (poly1, poly2, ctx);
    else
    {
        slong i;

        TEMPLATE(T, poly_fit_length) (poly1, len, ctx);

        for (i = 0; i < len; i++)
            TEMPLATE(T, set) (poly1->coeffs + i, poly2->coeffs + i, ctx);
        _TEMPLATE(T, poly_set_length) (poly1, len, ctx);
        _TEMPLATE(T, poly_normalise) (poly1, ctx);
    }
}


#endif
