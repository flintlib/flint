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

void elem_divexact_parent(elem_ptr res, elem_srcptr op1, elem_srcptr cont, const ring_t ring, const ring_t parent);

void
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


void elem_gcd_parent(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring, const ring_t parent);

void
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
        switch (ring->type)
        {
            case TYPE_FMPZ:
                fmpz_gcd(res, op1, op2);
                break;

            /* assumed to be a field */
            case TYPE_MOD:
                elem_one(res, ring);
                break;

            default:
                NOT_IMPLEMENTED("gcd_parent", ring);
        }
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

