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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

static __inline__ void
elem_poly_div_content_recursive(elem_poly_struct * poly, elem_srcptr cont,
        const ring_t cont_ring, const ring_t poly_ring)
{
    long i, size = RING_PARENT(poly_ring)->size;

    for (i = 0; i < poly->length; i++)
    {
        elem_div_content_recursive(INDEX(poly->coeffs, i, size), cont, cont_ring, RING_PARENT(poly_ring));
    }
}

void
elem_div_content_recursive(elem_ptr obj, elem_srcptr cont, const ring_t cont_ring, const ring_t obj_ring)
{
    if (cont_ring == obj_ring)
    {
        elem_divexact(obj, obj, cont, cont_ring);
    }
    else if (obj_ring->type == TYPE_POLY)
    {
        elem_poly_div_content_recursive(obj, cont, cont_ring, obj_ring);
    }
    else
    {
        NOT_IMPLEMENTED("div_content_recursive", obj_ring);
    }
}
