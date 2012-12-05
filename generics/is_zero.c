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
elem_is_zero(elem_srcptr x, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            return fmpz_is_zero(x);

        case TYPE_LIMB:
            return *((mp_srcptr) x) == 0;

        case TYPE_POLY:
            return ((elem_poly_struct *) x)->length == 0;

        case TYPE_MOD:
            return elem_is_zero(x, ring->parent);

        case TYPE_FRAC:
            return elem_is_zero(NUMER(x, ring), ring->numer);

        case TYPE_COMPLEX:
            return elem_is_zero(REALPART(x, ring), ring->parent) &&
                   elem_is_zero(IMAGPART(x, ring), ring->parent);

        default:
            NOT_IMPLEMENTED("is_zero", ring);
    }
}

