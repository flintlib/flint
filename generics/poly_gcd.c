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

int elem_leading_sign(elem_srcptr x, const ring_t ring);

#define ELEM_VEC_NORM(vec, len, ring) \
    do { \
        while ((len) > 0 && elem_is_zero(INDEX((vec), (len) - 1, (ring)->size), (ring))) \
            (len)--; \
    } while (0)


void
_elem_vec_gcd(elem_ptr res, elem_srcptr vec, long len, elem_srcptr input, const ring_t ring)
{
    long i, size = ring->size;

    if (len == 0)
    {
        elem_set(res, input, ring);
    }
    else
    {
        elem_gcd(res, vec, input, ring);

        for (i = 1; i < len && !elem_is_one(res, ring); i++)
            elem_gcd(res, res, INDEX(vec, i, size), ring);
    }
}

void
_elem_poly_pseudo_rem_cohen(elem_ptr R, elem_srcptr A, long lenA, elem_srcptr B, long lenB, const ring_t ring)
{
    long e, size = ring->size;
    elem_srcptr leadB = INDEX(B, lenB - 1, size);
    elem_ptr pow;

    if (lenB == 1)
    {
        _elem_vec_zero(R, lenA, ring);
        return;
    }

    if (R != A)
        _elem_vec_set(R, A, lenA, ring);

    e = lenA - lenB + 1;

    while (lenA >= lenB)
    {
        _elem_vec_scalar_mul(R, R, lenA - 1, leadB, ring);
        _elem_vec_scalar_submul(INDEX(R, lenA - lenB, size), B, lenB - 1, SRC_INDEX(R, lenA - 1, size), ring);
        elem_zero(INDEX(R, lenA - 1, size), ring);

        lenA--;
        ELEM_VEC_NORM(R, lenA, ring);

        e--;
    }

    ELEM_TMP_INIT(pow, ring);
    elem_pow_ui(pow, leadB, e, ring);
    _elem_vec_scalar_mul(R, R, lenA, pow, ring);
    ELEM_TMP_CLEAR(pow, ring);
}


void
_elem_poly_gcd_subresultant(elem_ptr res, elem_srcptr poly1, long len1, elem_srcptr poly2, long len2, const ring_t ring)
{
    if (len2 == 1)
    {
        _elem_vec_gcd(res, poly1, len1, poly2, ring);
    }
    else
    {
        elem_ptr a, b, d, g, h;
        elem_ptr A, B, W;
        long alloc, lenA, lenB, size = ring->size;

        alloc = len1 + len2 + 5;
        W = _elem_vec_init(alloc, ring);
        A = INDEX(W, 0, size);
        B = INDEX(A, len1, size);
        a = INDEX(B, len2, size);
        b = INDEX(a, 1, size);
        d = INDEX(b, 1, size);
        g = INDEX(d, 1, size);
        h = INDEX(g, 1, size);

        lenA = len1;
        lenB = len2;
        _elem_vec_gcd(a, poly1, lenA, a, ring);
        _elem_vec_gcd(b, poly2, lenB, b, ring);
        _elem_vec_scalar_divexact(A, poly1, lenA, a, ring);
        _elem_vec_scalar_divexact(B, poly2, lenB, b, ring);

        elem_gcd(d, a, b, ring);
        elem_one(g, ring);
        elem_one(h, ring);

        while (1)
        {
            const long delta = lenA - lenB;

            _elem_poly_pseudo_rem_cohen(A, A, lenA, B, lenB, ring);

            ELEM_VEC_NORM(A, lenA, ring);

            if (lenA <= 1)
                break;

            { /* Swap A and B */
                elem_ptr T;
                long len;
                T = A, A = B, B = T, len = lenA, lenA = lenB, lenB = len;
            }

            if (delta == 1)
            {
                elem_mul(b, g, h, ring);
                _elem_vec_scalar_divexact(B, B, lenB, b, ring);
                elem_set(g, SRC_INDEX(A, lenA - 1, size), ring);
                elem_set(h, g, ring);
            }
            else
            {
                elem_pow_ui(a, h, delta, ring);
                elem_mul(b, g, a, ring);
                _elem_vec_scalar_divexact(B, B, lenB, b, ring);
                elem_pow_ui(b, SRC_INDEX(A, lenA - 1, size), delta, ring);
                elem_mul(g, h, b, ring);
                elem_divexact(h, g, a, ring);
                elem_set(g, SRC_INDEX(A, lenA - 1, size), ring);
            }
        }

        if (lenA == 1)
        {
            elem_set(res, d, ring);
            _elem_vec_zero(INDEX(res, 1, size), len2 - 1, ring);
        }
        else
        {
            elem_zero(b, ring);
            _elem_vec_gcd(b, B, lenB, b, ring);
            _elem_vec_scalar_divexact(B, B, lenB, b, ring);

            if (elem_leading_sign(INDEX(B, lenB - 1, size), ring) < 0)
                elem_neg(d, d, ring);

            _elem_vec_scalar_mul(res, B, lenB, d, ring);

            if (len2 >= lenB)
                _elem_vec_zero(INDEX(res, lenB, size), len2 - lenB, ring);
        }

        _elem_vec_clear(W, alloc, ring);
    }
}

void
elem_poly_gcd_subresultant(elem_poly_t res, const elem_poly_t poly1, const elem_poly_t poly2, const ring_t ring)
{
    if (poly1->length < poly2->length)
    {
        elem_poly_gcd_subresultant(res, poly2, poly1, ring);
    }
    else /* len1 >= len2 >= 0 */
    {
        const long len1 = poly1->length;
        const long len2 = poly2->length;
        
        if (len1 == 0) /* len1 = len2 = 0 */
        {
            elem_zero(res, ring);
        }
        else if (len2 == 0) /* len1 > len2 = 0 */
        {
            if (elem_leading_sign(poly1, ring) > 0)
                elem_poly_set(res, poly1, ring);
            else
                elem_poly_neg(res, poly1, ring);
        }
        else /* len1 >= len2 >= 1 */
        {
            /* underscore code automatically handles aliasing */
           
            elem_poly_fit_length(res, len2, ring);
            _elem_poly_gcd_subresultant(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, RING_PARENT(ring));
            elem_poly_set_length(res, len2, ring);
            elem_poly_normalise(res, ring);
        }
    }
}
