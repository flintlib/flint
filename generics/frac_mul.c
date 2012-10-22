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
elem_frac_mul(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring)
{
    if (elem_is_zero(op1, ring) || elem_is_zero(op2, ring))
    {
        elem_zero(res, ring);
    }
    else if (res == op1 || res == op2)
    {
        elem_ptr t;
        ELEM_TMP_INIT(t, ring);
        elem_frac_mul(t, op1, op2, ring);
        elem_swap(res, t, ring);
        ELEM_TMP_CLEAR(t, ring);
    }
    else
    {
        elem_mul(NUMER(res, ring), NUMER(op1, ring), NUMER(op2, ring), RING_NUMER(ring));
        elem_mul(DENOM(res, ring), DENOM(op1, ring), DENOM(op2, ring), RING_DENOM(ring));
        elem_frac_canonicalise(res, ring);
    }
}

