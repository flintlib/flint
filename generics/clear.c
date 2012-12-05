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
elem_poly_clear(elem_poly_struct * poly, const ring_t ring)
{
    long i, size;
    size = RING_PARENT(ring)->size;

    for (i = 0; i < poly->alloc; i++)
        elem_clear(INDEX(poly->coeffs, i, size), ring->parent);

    flint_free(poly->coeffs);
}

void
elem_clear(elem_ptr elem, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_clear((fmpz *) elem);
            break;

        case TYPE_LIMB:
            break;

        case TYPE_POLY:
            elem_poly_clear((elem_poly_struct *) elem, ring);
            break;

        case TYPE_MOD:
            elem_clear(elem, ring->parent);
            break;

        case TYPE_FRAC:
            elem_clear(NUMER(elem, ring), ring->numer);
            elem_clear(DENOM(elem, ring), ring->denom);
            break;

        case TYPE_COMPLEX:
            elem_clear(REALPART(elem, ring), ring->parent);
            elem_clear(IMAGPART(elem, ring), ring->parent);
            break;

        default:
            NOT_IMPLEMENTED("clear", ring);
    }
}

