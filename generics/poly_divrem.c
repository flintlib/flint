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

static int
euclidean_test(elem_srcptr X, elem_srcptr L, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            return fmpz_cmpabs(X, L) < 0;

        case TYPE_MOD:
            return elem_is_zero(X, ring);

        case TYPE_POLY:
            return ((const elem_poly_struct *) X)->length < ((const elem_poly_struct *) L)->length;

        case TYPE_FRAC:
            /* assume that it's a field */
            if (RING_NUMER(ring) == RING_DENOM(ring))
                return elem_is_zero(X, ring);
            else
                return euclidean_test(NUMER(X, ring), NUMER(L, ring), RING_NUMER(ring));

        case TYPE_LIMB:
        default:
            NOT_IMPLEMENTED("euclidean test", ring); 
    }
}



void _elem_poly_divrem(elem_ptr Q, elem_ptr R, elem_srcptr A,
    long lenA, elem_srcptr B, long lenB, const ring_t ring)
{
    elem_ptr tmp;
    elem_srcptr leadB;
    long i, iQ, iR, size;

    if (ring->type == TYPE_FMPZ && 0)
    {
        _fmpz_poly_divrem(Q, R, A, lenA, B, lenB);
        return;
    }

    size = ring->size;
    leadB = SRC_INDEX(B, lenB - 1, size);
    ELEM_TMP_INIT(tmp, ring);

    if (R != A)
    {
        _elem_vec_set(R, A, lenA, ring);
    }

    for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
    {
        if (euclidean_test(INDEX(R, iR, size), leadB, ring))
        {
            elem_zero(INDEX(Q, iQ, size), ring);
        }
        else
        {
            /* TODO: this could use div, not writing out the remainder */
            elem_divrem(INDEX(Q, iQ, size), tmp, SRC_INDEX(R, iR, size), leadB, ring);

            for (i = 0; i < lenB; i++)
            {
                /* R[iQ + i] -= B[i] * Q[iQ] */
                elem_mul(tmp, INDEX(B, i, size), SRC_INDEX(Q, iQ, size), ring);
                elem_sub(INDEX(R, iQ + i, size), SRC_INDEX(R, iQ + i, size), tmp, ring);
            }
        }
    }

    ELEM_TMP_CLEAR(tmp, ring);
}

void
elem_poly_divrem(elem_poly_t Q, elem_poly_t R,
                      const elem_poly_t A, const elem_poly_t B, const ring_t ring)
{
    long lenA = A->length;
    long lenB = B->length;

    if (lenB == 0)
    {
        printf("Exception: division by zero in elem_poly_divrem\n");
        abort();
    }

    if (lenA < lenB)
    {
        elem_set(R, A, ring);
        elem_zero(Q, ring);
        return;
    }

    if (Q == A || Q == B)
    {
        elem_poly_t tmp;
        elem_init(tmp, ring);
        elem_poly_divrem(tmp, R, A, B, ring);
        elem_swap(Q, tmp, ring);
        elem_clear(tmp, ring);
        return;
    }

    if (R == B)
    {
        elem_poly_t tmp;
        elem_init(tmp, ring);
        elem_poly_divrem(Q, tmp, A, B, ring);
        elem_swap(R, tmp, ring);
        elem_clear(tmp, ring);
        return;
    }

    /* TODO: lenB - 1, if a field... */
    elem_poly_fit_length(R, lenA, ring);
    elem_poly_fit_length(Q, lenA - lenB + 1, ring);

    _elem_poly_divrem(Q->coeffs, R->coeffs, A->coeffs, lenA, B->coeffs, lenB, ring->parent);

    elem_poly_set_length(Q, lenA - lenB + 1, ring);
    elem_poly_set_length(R, lenA, ring);
    elem_poly_normalise(Q, ring);
    elem_poly_normalise(R, ring);
}

