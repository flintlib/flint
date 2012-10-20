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
elem_poly_init(elem_poly_struct * poly, const ring_t ring)
{
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void
elem_init(elem_ptr elem, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_init(elem);
            break;

        case TYPE_LIMB:
            *((mp_ptr) elem) = 0;
            break;

        case TYPE_POLY:
            elem_poly_init(elem, ring);
            break;

        case TYPE_MOD:
            elem_init(elem, ring->parent);
            break;

        default:
            NOT_IMPLEMENTED("init", ring);
    }
}

void
gen_init(gen_t x, const ring_t ring)
{
    x->ring = (ring_struct *) ring;
    x->elem = flint_malloc(ring->size);
    elem_init(x->elem, ring);
}
