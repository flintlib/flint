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
    Copyright (C) 2013 William Hart

******************************************************************************/

#include "generics.h"

static __inline__ void
elem_poly_content_recursive(elem_ptr cont, const elem_poly_struct * poly,
        const ring_t cont_ring, const ring_t poly_ring)
{
    long i, size = RING_PARENT(poly_ring)->size;

    /* todo: should break when == 1 */
    for (i = 0; i < poly->length; i++)
    {
        elem_content_recursive(cont, INDEX(poly->coeffs, i, size), cont_ring, RING_PARENT(poly_ring));
    }
}

void
elem_content_recursive(elem_ptr cont, elem_srcptr obj, const ring_t cont_ring, const ring_t obj_ring)
{
    if (cont_ring == obj_ring)
    {
        switch (cont_ring->type)
        {
            case TYPE_FMPZ:
                fmpz_gcd(cont, cont, obj);
                break;

            case TYPE_MPZ:
                mpz_gcd(cont, cont, obj);
                break;

            /* assumed to be a field */
            case TYPE_MOD:
                elem_set_si(cont, 1, obj_ring);
/*                elem_one(cont, obj_ring); */
                break;

            default:
                NOT_IMPLEMENTED("content_recursive", obj_ring);
        }
    }
    else
    {
        switch (obj_ring->type)
        {
            case TYPE_POLY:
                elem_poly_content_recursive(cont, obj, cont_ring, obj_ring);
                break;

            /* assumed to be a field */
            default:
                NOT_IMPLEMENTED("content_recursive", obj_ring);
        }
    }
}
