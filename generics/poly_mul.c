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
_elem_poly_mul(elem_ptr res, elem_srcptr poly1, long len1,
    elem_srcptr poly2, long len2, const ring_t ring)
{
/*
    if (ring->type == TYPE_MOD && RING_PARENT(ring)->type == TYPE_FMPZ)
    {
        _fmpz_mod_poly_mul(res, poly1, len1, poly2, len2, &RING_MODULUS(ring)->z);
        return;
    }
*/
    long size = ring->size;

    if (len1 == 1 && len2 == 1)
    {
        elem_mul(res, poly1, poly2, ring);
    }
    else
    {
        long i;

        _elem_vec_scalar_mul(res, poly1, len1, poly2, ring);

        _elem_vec_scalar_mul(INDEX(res, len1, size),
                            SRC_INDEX(poly2, 1, size), len2 - 1,
                            SRC_INDEX(poly1, len1 - 1, size), ring);

        for (i = 0; i < len1 - 1; i++)
            _elem_vec_scalar_addmul(INDEX(res, i + 1, size),
                    SRC_INDEX(poly2, 1, size), len2 - 1,
                    SRC_INDEX(poly1, i, size), ring);
    }
}
