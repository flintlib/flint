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

void
elem_zero(elem_ptr x, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_zero(x);
            break;

        case TYPE_LIMB:
            *((mp_ptr) x) = 0;
            break;

        case TYPE_POLY:
            ((elem_poly_struct *) x)->length = 0;
            break;

        case TYPE_MOD:
            elem_zero(x, ring->parent);
            break;

        case TYPE_FRAC:
            elem_zero(NUMER(x, ring), RING_NUMER(ring));
            elem_one(DENOM(x, ring), RING_DENOM(ring));
            break;

        default:
            NOT_IMPLEMENTED("zero", ring);
    }
}

void
gen_zero(gen_t x)
{
    elem_zero(x->elem, x->ring);
}
