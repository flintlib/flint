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
elem_set(elem_ptr res, elem_srcptr src, const ring_t ring)
{
    if (res != src)
    {
        switch (ring->type)
        {
            case TYPE_FMPZ:
                fmpz_set(res, src);
                break;

            case TYPE_LIMB:
                *((mp_ptr) res) = *((mp_srcptr) src);
                break;

            case TYPE_MOD:
                elem_set(res, src, ring->parent);
                break;

            case TYPE_POLY:
                elem_poly_set(res, src, ring);
                break;

            case TYPE_FRAC:
                elem_set(NUMER(res, ring), NUMER(src, ring), RING_NUMER(ring));
                elem_set(DENOM(res, ring), DENOM(src, ring), RING_DENOM(ring));
                break;

            case TYPE_COMPLEX:
                elem_set(REALPART(res, ring), REALPART(src, ring), RING_PARENT(ring));
                elem_set(IMAGPART(res, ring), IMAGPART(src, ring), RING_PARENT(ring));
                break;

            default:
                NOT_IMPLEMENTED("set", ring);
        }
    }
}

