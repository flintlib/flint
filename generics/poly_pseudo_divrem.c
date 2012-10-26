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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

void
_elem_poly_pseudo_divrem(elem_ptr Q, elem_ptr R, ulong * d,
    elem_srcptr A, long lenA, elem_srcptr B, long lenB, const ring_t ring)
{
    long iQ = lenA - lenB, iR = lenA - 1;
    long size = ring->size;
    elem_srcptr leadB = SRC_INDEX(B, lenB - 1, size);
    elem_ptr rem;

    ELEM_TMP_INIT(rem, ring);

    *d = 0;
    _elem_vec_zero(Q, lenA - lenB + 1, ring);

    if (R != A)
        _elem_vec_set(R, A, lenA, ring);

    while (iR >= lenB - 1)
    {
        elem_divrem(INDEX(Q, iQ, size), rem, SRC_INDEX(R, iR, size), leadB, ring);

        if (!elem_is_zero(rem, ring))
        {
            _elem_vec_scalar_mul(Q, Q, lenA - lenB + 1, leadB, ring);
            elem_set(INDEX(Q, iQ, size), INDEX(R, iR, size), ring);
            _elem_vec_scalar_mul(R, R, lenA, leadB, ring);
            (*d)++;
        }

        if (lenB > 1)
            _elem_vec_scalar_submul(INDEX(R, iR - lenB + 1, size),
                B, lenB - 1, SRC_INDEX(Q, iQ, size), ring);

        elem_zero(INDEX(R, iR, size), ring);

        iR--;
        iQ--;
    }

    ELEM_TMP_CLEAR(rem, ring);
}


void
elem_poly_pseudo_divrem(elem_poly_t Q, elem_poly_t R,
                                 ulong * d, const elem_poly_t A,
                                 const elem_poly_t B, const ring_t ring)
{
    long lenq, lenr;

    if (B->length == 0)
    {
        printf("Exception: division by zero in elem_poly_pseudo_divrem\n");
        abort();
    }
    else if (Q == R)
    {
        printf("Exception: output arguments Q and R may not be aliased\n");
        abort();
    }
    else if (A->length < B->length)
    {
        elem_set(R, A, ring);
        elem_zero(Q, ring);
        *d = 0;
    }
    else if (Q == A || Q == B)
    {
        elem_poly_t tmp;
        elem_init(tmp, ring);
        elem_poly_pseudo_divrem(tmp, R, d, A, B, ring);
        elem_swap(Q, tmp, ring);
        elem_clear(tmp, ring);
    }
    else if (R == B)
    {
        elem_poly_t tmp;
        elem_init(tmp, ring);
        elem_poly_pseudo_divrem(Q, tmp, d, A, B, ring);
        elem_swap(R, tmp, ring);
        elem_clear(tmp, ring);
    }
    else
    {
        lenq = A->length - B->length + 1;
        lenr = A->length;

        elem_poly_fit_length(Q, lenq, ring);
        elem_poly_fit_length(R, lenr, ring);

        _elem_poly_pseudo_divrem(Q->coeffs, R->coeffs, d, A->coeffs, A->length,
            B->coeffs, B->length, RING_PARENT(ring));

        elem_poly_set_length(Q, lenq, ring);
        elem_poly_set_length(R, lenr - 1, ring);
        elem_poly_normalise(R, ring);
    }
}

void
gen_pseudo_divrem(gen_t Q, gen_t R, ulong * d, const gen_t A, const gen_t B)
{
    if (Q->ring == A->ring && A->ring == B->ring && R->ring == A->ring && Q->ring->type == TYPE_POLY)
    {
        elem_poly_pseudo_divrem(Q->elem, R->elem, d, A->elem, B->elem, Q->ring);
    }
    else
    {
        NOT_IMPLEMENTED("gen_pseudo_divrem coercing into ", Q->ring);
    }
}

