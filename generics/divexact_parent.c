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

static void
elem_poly_divexact_parent(elem_poly_struct * res,
    const elem_poly_struct * op1, elem_srcptr cont,
    const ring_t ring, const ring_t parent)
{
    long i, len = op1->length, size = RING_PARENT(ring)->size;

    _elem_poly_fit_length(res, len, ring);

    for (i = 0; i < len; i++)
    {
        elem_divexact_parent(INDEX(res->coeffs, i, size),
            SRC_INDEX(op1->coeffs, i, size),
            cont, RING_PARENT(ring), parent);
    }

    _elem_poly_set_length(res, len, ring);
}

void
elem_divexact_parent(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring, const ring_t parent)
{
    if (ring == parent)
    {
        elem_divexact(res, op1, op2, ring);
    }
    else if (ring->type == TYPE_POLY)
    {
        elem_poly_divexact_parent(res, op1, op2, ring, parent);
    }
    else
    {
        NOT_IMPLEMENTED("divexact_parent", ring);
    }
}

