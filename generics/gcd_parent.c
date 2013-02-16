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

static void
elem_poly_gcd_parent(elem_ptr res, const elem_poly_struct * op1, elem_srcptr op2, const ring_t ring, const ring_t parent)
{
    long i, len = op1->length, size = RING_PARENT(ring)->size;

    if (len == 0)
    {
        elem_set(res, op2, parent);
    }
    else
    {
        elem_gcd_parent(res, op1->coeffs, op2, RING_PARENT(ring), parent);

        for (i = 1; i < len && !elem_is_one(res, parent); i++)
            elem_gcd_parent(res, INDEX(op1->coeffs, i, size), res, RING_PARENT(ring), parent);
    }
}

void
elem_gcd_parent(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring, const ring_t parent)
{
    if (ring == parent)
    {
        elem_gcd(res, op1, op2, ring);
    }
    else
    {
        switch (ring->type)
        {
            case TYPE_POLY:
                elem_poly_gcd_parent(res, op1, op2, ring, parent);
                break;

            /* assumed to be a field */
            default:
                NOT_IMPLEMENTED("gcd_parent (recursive)", ring);
        }

    }
}

void
elem_gcd(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_gcd(res, op1, op2);
            break;

        case TYPE_MPZ:
            mpz_gcd(res, op1, op2);
            break;

        /* assumed to be a field */
        case TYPE_MOD:
            elem_one(res, ring);
            break;

        case TYPE_POLY:
            elem_poly_gcd_subresultant(res, op1, op2, ring);
            break;

        default:
            NOT_IMPLEMENTED("gcd", ring);
    }
}
