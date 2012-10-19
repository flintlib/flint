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
_elem_poly_print(elem_srcptr poly, long len, const ring_t ring)
{
    long i;

    printf("[");

    for (i = 0; i < len; i++)
    {
        elem_print(poly + i, ring);

        if (i < len - 1)
            printf(", ");
    }

    printf("]");
}

void
elem_print(const elem_t elem, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_print(&elem->z);
            break;

        case TYPE_LIMB:
            printf("%lu", elem->n);
            break;

        case TYPE_POLY:
            _elem_poly_print(elem->poly->coeffs, elem->poly->length, ring->parent);
            break;

        case TYPE_MOD:
            elem_print(elem, ring->parent);
            break;

        default:
            NOT_IMPLEMENTED("print", ring);
    }
}

void gen_print(gen_t x)
{
    printf("element of ");
    ring_print(x->ring);
    printf(":\n");
    elem_print(&x->elem, x->ring);
    printf("\n");
}
