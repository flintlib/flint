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

/*
divrem over (R[x] / R)
*/

void
_elem_frac_poly_divrem(elem_ptr Q, elem_ptr Qden,
    elem_ptr R, elem_ptr Rden,
    elem_srcptr A, elem_srcptr Aden, long lenA,
    elem_srcptr B, elem_srcptr Bden, long lenB, const ring_t ring)
{
    long lenQ = lenA - lenB + 1;
    long lenR = lenB - 1;
    long size = ring->size;
    elem_ptr den;
    elem_srcptr leadB = SRC_INDEX(B, lenB - 1, size);
    ulong d;

    /* TODO: lenB == 1 case efficiently */
    _elem_poly_pseudo_divrem(Q, R, &d, A, lenA, B, lenB, ring);

    for ( ; lenR != 0 &&
        elem_is_zero(SRC_INDEX(R, lenR - 1, size), ring); lenR--);

    /* TODO: lead^d = +/- 1 and other cases more efficiently */
    /* see also fmpq_poly_divrem, which avoids the post-canonicalisation */

    ELEM_TMP_INIT(den, ring);

    elem_pow_ui(den, leadB, d, ring);
    elem_mul(den, den, Aden, ring);
    elem_set(Rden, den, ring);

    elem_set(Qden, den, ring);

    _elem_vec_scalar_mul(Q, Q, lenQ, Bden, ring);

    ELEM_TMP_CLEAR(den, ring);
}

void
_elem_frac_poly_divrem_wrap(elem_poly_struct * Q, elem_ptr Qden,
    elem_poly_struct * R, elem_ptr Rden,
    const elem_poly_struct * A, elem_srcptr Aden,
    const elem_poly_struct * B, elem_srcptr Bden, const ring_t num_ring, const ring_t den_ring)
{
    long lenA, lenB, lenQ, lenR;

    lenA = A->length;
    lenB = B->length;
    lenQ = lenA - lenB + 1;
    lenR = lenB - 1;

    elem_poly_fit_length(Q, lenQ, num_ring);
    elem_poly_fit_length(R, lenA, num_ring); /* XXX: needed space for pseudo division */

    _elem_frac_poly_divrem(Q->coeffs, Qden, R->coeffs, Rden,
        A->coeffs, Aden, A->length,
        B->coeffs, Bden, B->length, RING_PARENT(num_ring));

    elem_poly_set_length(Q, lenQ, num_ring);
    elem_poly_set_length(R, lenR, num_ring);
    elem_poly_normalise(R, num_ring);

    _elem_frac_canonicalise(Q, Qden, num_ring, den_ring);
    _elem_frac_canonicalise(R, Rden, num_ring, den_ring);
}

void
elem_frac_poly_divrem(elem_ptr Q, elem_ptr R, elem_srcptr A, elem_srcptr B, const ring_t ring)
{
    if (elem_is_zero(B, ring))
    {
        printf("Exception: division by zero in elem_frac_poly_divrem\n");
        abort();
    }
    else if (Q == R)
    {
        printf("Exception: output arguments Q and R may not be aliased\n");
        abort();
    }
    else if (((const elem_poly_struct *) A)->length < ((const elem_poly_struct *) B)->length)
    {
        elem_set(R, A, ring);
        elem_zero(Q, ring);
    }
    else if (Q == A || Q == B)
    {
        elem_ptr tmp;
        ELEM_TMP_INIT(tmp, ring);
        elem_frac_poly_divrem(tmp, R, A, B, ring);
        elem_swap(Q, tmp, ring);
        ELEM_TMP_CLEAR(tmp, ring);
    }
    else if (R == B)
    {
        elem_ptr tmp;
        ELEM_TMP_INIT(tmp, ring);
        elem_frac_poly_divrem(Q, tmp, A, B, ring);
        elem_swap(R, tmp, ring);
        ELEM_TMP_CLEAR(tmp, ring);
    }
    else
    {
        /* if the denominators are the same as coefficients, we can do pseudo division */
        if (RING_DENOM(ring) == RING_PARENT(RING_NUMER(ring)))
        {
            _elem_frac_poly_divrem_wrap(NUMER(Q, ring), DENOM(Q, ring),
                NUMER(R, ring), DENOM(R, ring),
                NUMER(A, ring), DENOM(A, ring),
                NUMER(B, ring), DENOM(B, ring), RING_NUMER(ring), RING_DENOM(ring));
        }
        else
        {
            NOT_IMPLEMENTED("polynomial divrem with different denominator type", ring);
        }
    }
}

