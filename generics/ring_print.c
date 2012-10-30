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
ring_print(const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            printf("fmpzs");
            break;

        case TYPE_LIMB:
            printf("limbs");
            break;

        case TYPE_POLY:
            printf("polynomials over (");
            ring_print(ring->parent);
            printf(")");
            break;

        case TYPE_MOD:
            ring_print(ring->parent);
            printf(" modulo ");
            elem_print(RING_MODULUS(ring), ring->parent);
            break;

        case TYPE_FRAC:
            printf("(");
            ring_print(ring->numer);
            printf(") divided by (");
            ring_print(ring->denom);
            printf(")");
            break;

        case TYPE_MAT:
            printf("matrices over (");
            ring_print(ring->parent);
            printf(")");
            break;

        case TYPE_COMPLEX:
            printf("complex extension of (");
            ring_print(ring->parent);
            printf(")");
            break;

        default:
            printf("unknown ring\n");
            abort();
    }
}
