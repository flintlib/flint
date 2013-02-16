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

void
elem_swap(elem_ptr res, elem_ptr src, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_swap(res, src);
            break;

        case TYPE_MPZ:
            mpz_swap(res, src);
            break;

        case TYPE_LIMB:
            {
                mp_limb_t t = *((mp_ptr) res);
                *((mp_ptr) res) = *((mp_srcptr) src);
                *((mp_ptr) src) = t;
            }
            break;

        case TYPE_MOD:
            elem_swap(res, src, ring->parent);
            break;

        case TYPE_POLY:
            elem_poly_swap(res, src);
            break;

        case TYPE_FRAC:
            elem_swap(NUMER(res, ring), NUMER(src, ring), RING_NUMER(ring));
            elem_swap(DENOM(res, ring), DENOM(src, ring), RING_DENOM(ring));
            break;

        case TYPE_COMPLEX:
            elem_swap(REALPART(res, ring), REALPART(src, ring), RING_PARENT(ring));
            elem_swap(IMAGPART(res, ring), IMAGPART(src, ring), RING_PARENT(ring));
            break;

        default:
            NOT_IMPLEMENTED("swap", ring);
    }
}

