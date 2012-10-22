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

void elem_mul_parent(elem_ptr res, elem_srcptr op1, elem_srcptr cont, const ring_t ring, const ring_t parent);

/* todo: handle zero, efficient in-place...? */
void
elem_poly_mul_parent(elem_poly_struct * res, const elem_poly_struct * op1, elem_srcptr op2, const ring_t ring, const ring_t parent)
{
    long i, size = RING_PARENT(ring)->size;

    _elem_poly_fit_length(res, op1->length, ring);

    for (i = 0; i < op1->length; i++)
    {
        elem_mul_parent(INDEX(res->coeffs, i, size), SRC_INDEX(op1->coeffs, i, size), op2, RING_PARENT(ring), parent);
    }

    _elem_poly_set_length(res, op1->length, ring);
}

void
elem_mul_parent(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring, const ring_t parent)
{
    if (ring == parent)
    {
        elem_mul(res, op1, op2, ring);
    }
    else if (ring->type == TYPE_POLY)
    {
        elem_poly_mul_parent(res, op1, op2, ring, parent);
    }
    else
    {
        NOT_IMPLEMENTED("mul_parent", ring);
    }
}

void
elem_frac_add(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring)
{
    if (elem_is_zero(op1, ring))
    {
        elem_set(res, op2, ring);
    }
    else if (elem_is_zero(op2, ring))
    {
        elem_set(res, op1, ring);
    }
    else
    {
        elem_ptr t, u;

        ELEM_TMP_INIT(t, RING_NUMER(ring));
        ELEM_TMP_INIT(u, RING_NUMER(ring));

        elem_mul_parent(t, NUMER(op1, ring), DENOM(op2, ring), RING_NUMER(ring), RING_DENOM(ring));
        elem_mul_parent(u, NUMER(op2, ring), DENOM(op1, ring), RING_NUMER(ring), RING_DENOM(ring));

        elem_add(NUMER(res, ring), t, u, RING_NUMER(ring));
        elem_mul(DENOM(res, ring), DENOM(op1, ring), DENOM(op2, ring), RING_DENOM(ring));

        elem_frac_canonicalise(res, ring);

        ELEM_TMP_CLEAR(t, RING_NUMER(ring));
        ELEM_TMP_CLEAR(u, RING_NUMER(ring));
    }
}

