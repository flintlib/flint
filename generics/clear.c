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
elem_clear(elem_t elem, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_clear(&elem->z);
            break;

        case TYPE_LIMB:
            break;

        case TYPE_POLY:
            {
                long i;
                for (i = 0; i < elem->poly->alloc; i++)
                    elem_clear(((elem_ptr) elem->poly->coeffs) + i, ring->parent);
                flint_free(elem->poly);
            }
            break;

        case TYPE_MOD:
            elem_clear(elem, ring->parent);
            break;

        default:
            NOT_IMPLEMENTED("clear", ring);
    }
}

void
gen_clear(gen_t x)
{
    elem_clear(&x->elem, x->ring);
}
