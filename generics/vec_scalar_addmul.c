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
_elem_vec_scalar_addmul(elem_ptr res, elem_srcptr vec, long len, elem_srcptr c, const ring_t ring)
{
    long i, size = ring->size;
    elem_ptr t;

    /* XXX: should use tmp alloc */
    t = flint_malloc(ring->size);
    elem_init(t, ring);

    for (i = 0; i < len; i++)
    {
        elem_mul(t, SRC_INDEX(vec, i, size), c, ring);
        elem_add(INDEX(res, i, size), SRC_INDEX(res, i, size), t, ring);
    }

    elem_clear(t, ring);
    flint_free(t);
}
