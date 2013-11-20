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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
_TEMPLATE3(T, vec_scalar_submul, T) (TEMPLATE(T, struct) * poly1,
                                     const TEMPLATE(T, struct) * poly2,
                                     slong len2, const TEMPLATE(T, t) x,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    TEMPLATE(T, t) y;

    TEMPLATE(T, init) (y, ctx);

    for (i = 0; i < len2; i++)
    {
        TEMPLATE(T, mul) (y, poly2 + i, x, ctx);
        TEMPLATE(T, sub) (poly1 + i, poly1 + i, y, ctx);
    }

    TEMPLATE(T, clear) (y, ctx);
}


#endif
