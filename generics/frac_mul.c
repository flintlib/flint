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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

/* assumes: nonzero, no aliasing */
void
_elem_frac_mul(elem_ptr ropnum, elem_ptr ropden,
        elem_srcptr op1num, elem_srcptr op1den,
        elem_srcptr op2num, elem_srcptr op2den,
        const ring_t num_ring, const ring_t den_ring)
{
    elem_ptr t, u, v, a, b;

    /* likely case (e.g. both integral) */
    if (elem_equal(op1den, op2den, den_ring))
    {
        elem_mul(ropnum, op1num, op2num, num_ring);
        elem_mul(ropden, op1den, op2den, den_ring);
        return;
    }

    /* since there is no aliasing, we can use ropden as a temporary */
    a = ropden;
    elem_gcd_parent(a, op1num, op2den, num_ring, den_ring);

    if (elem_is_one(a, den_ring))
    {
        elem_gcd_parent(a, op2num, op1den, num_ring, den_ring);

        if (elem_is_one(a, den_ring))
        {
            elem_mul(ropnum, op1num, op2num, num_ring);
            elem_mul(ropden, op1den, op2den, den_ring);
        }
        else
        {
            elem_divexact_parent(ropnum, op2num, a, num_ring, den_ring);
            elem_mul(ropnum, op1num, ropnum, num_ring);
            elem_divexact(ropden, op1den, a, den_ring);
            elem_mul(ropden, ropden, op2den, den_ring);
        }
    }
    else
    {
        ELEM_TMP_INIT(b, den_ring);

        elem_gcd_parent(b, op2num, op1den, num_ring, den_ring);

        if (elem_is_one(b, den_ring))
        {
            elem_divexact_parent(ropnum, op1num, a, num_ring, den_ring);
            elem_mul(ropnum, ropnum, op2num, num_ring);
            elem_divexact(ropden, op2den, a, den_ring);
            elem_mul(ropden, op1den, ropden, den_ring);
        }
        else
        {
            ELEM_TMP_INIT(t, num_ring);
            ELEM_TMP_INIT(u, num_ring);

            elem_divexact_parent(t, op1num, a, num_ring, den_ring);
            elem_divexact_parent(u, op2num, b, num_ring, den_ring);
            elem_mul(ropnum, t, u, num_ring);

            ELEM_TMP_CLEAR(t, num_ring);
            ELEM_TMP_CLEAR(u, num_ring);

            ELEM_TMP_INIT(v, den_ring);

            elem_divexact(v, op1den, b, den_ring);
            elem_divexact(b, op2den, a, den_ring);
            elem_mul(ropden, v, b, den_ring);

            ELEM_TMP_CLEAR(v, den_ring);
        }

        ELEM_TMP_CLEAR(b, den_ring);
    }
}

void
elem_frac_mul(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring)
{
/*
    if (RING_NUMER(ring)->type == TYPE_FMPZ)
    {
        fmpq_mul(res, op1, op2);
        return;
    }
*/

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
        _elem_frac_mul(NUMER(res, ring), DENOM(res, ring),
            NUMER(op1, ring), DENOM(op1, ring),
            NUMER(op2, ring), DENOM(op2, ring),
            RING_NUMER(ring), RING_DENOM(ring));

/*
    naive version (but sometimes faster):
    elem_mul(NUMER(res, ring), NUMER(op1, ring), NUMER(op2, ring), RING_NUMER(ring));
    elem_mul(DENOM(res, ring), DENOM(op1, ring), DENOM(op2, ring), RING_DENOM(ring));
    elem_frac_canonicalise(res, ring);
*/
    }
}

