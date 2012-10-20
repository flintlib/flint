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
_elem_poly_fit_length(elem_poly_struct * poly, long len, const ring_t poly_ring)
{
    long i;

    if (len > poly->alloc)
    {
        elem_ptr ptr;
        long size = RING_PARENT(poly_ring)->size;

        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;

        ptr = poly->coeffs = flint_realloc(poly->coeffs, len * size);

        for (i = poly->alloc; i < len; i++)
            elem_init(INDEX(ptr, i, size), poly_ring->parent);

        poly->alloc = len;
    }
}
