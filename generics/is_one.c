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

int
elem_is_one(elem_srcptr x, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            return fmpz_is_one(x);

        case TYPE_LIMB:
            return *((mp_srcptr) x) == 1;

        /* XXX: mod 1? */
        case TYPE_MOD:
            return elem_is_one(x, ring->parent);

        case TYPE_FRAC:
            return elem_is_one(NUMER(x, ring), ring->numer) && 
                   elem_is_one(DENOM(x, ring), ring->denom);

        case TYPE_POLY:
            {
                const elem_poly_struct * poly = x;
                return (poly->length == 1) && elem_is_one(poly->coeffs, RING_PARENT(ring));
            }

        default:
            NOT_IMPLEMENTED("is_one", ring);
    }
}

int
gen_is_one(const gen_t x)
{
    return elem_is_one(x->elem, x->ring);
}
