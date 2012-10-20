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

/* XXX: rename? */
void
elem_set_coeff_si(elem_poly_struct * poly, long index, long value, const ring_t ring)
{
    long i, size = RING_PARENT(ring)->size;
    elem_ptr ptr;

    _elem_poly_fit_length(poly, index + 1, ring);
    ptr = poly->coeffs;

    for (i = poly->length; i < index; i++)
        elem_zero(INDEX(ptr, i, size), ring->parent);

    elem_set_si(INDEX(ptr, index, size), value, ring->parent);

    _elem_poly_set_length(poly, FLINT_MAX(poly->length, index + 1), ring);
    _elem_poly_normalise(poly, ring);
}

void
gen_set_coeff_si(gen_t x, long index, long value)
{
    if (x->ring->type != TYPE_POLY)
    {
        printf("set_coeff_si: ring must be a polynomial ring\n");
        abort();
    }

    elem_set_coeff_si(x->elem, index, value, x->ring);
}
