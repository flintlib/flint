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

int
elem_leading_sign(elem_srcptr x, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            return fmpz_sgn(x);

        case TYPE_POLY:
            {
                const elem_poly_struct * poly = x;
                if (poly->length == 0)
                    return 0;
                else
                    return elem_leading_sign(INDEX(poly->coeffs, poly->length - 1,
                        RING_PARENT(ring)->size), RING_PARENT(ring));
            }

        default:
            NOT_IMPLEMENTED("leading_sign", ring);
    }
}

void
elem_frac_canonicalise(elem_srcptr x, const ring_t ring)
{
    elem_ptr g;

    g = flint_malloc(RING_DENOM(ring)->size);
    elem_init(g, RING_DENOM(ring));

    elem_content_recursive(g, NUMER(x, ring), RING_DENOM(ring), RING_NUMER(ring));
    elem_content_recursive(g, DENOM(x, ring), RING_DENOM(ring), RING_DENOM(ring));

    if (!elem_is_one(g, RING_DENOM(ring)))
    {
        elem_div_content_recursive(NUMER(x, ring), g, RING_DENOM(ring), RING_NUMER(ring));
        elem_div_content_recursive(DENOM(x, ring), g, RING_DENOM(ring), RING_DENOM(ring));
    }

    if (elem_leading_sign(DENOM(x, ring), RING_DENOM(ring)) < 0)
    {
        elem_neg(NUMER(x, ring), NUMER(x, ring), RING_NUMER(ring));
        elem_neg(DENOM(x, ring), DENOM(x, ring), RING_DENOM(ring));
    }

    elem_clear(g, RING_DENOM(ring));
    flint_free(g);
}
