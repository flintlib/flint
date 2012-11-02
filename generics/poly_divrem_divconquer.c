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


#define FLINT_DIVREM_DIVCONQUER_CUTOFF 4

#define POS(vec, i, ring) INDEX(vec, i, (ring)->size)

void
_elem_poly_divrem_divconquer_recursive(elem_ptr Q, elem_ptr BQ, elem_ptr W, elem_srcptr A, elem_srcptr B, long lenB, const ring_t ring)
{
    if (lenB <= FLINT_DIVREM_DIVCONQUER_CUTOFF)
    {
        _elem_vec_zero(BQ, lenB - 1, ring);
        _elem_vec_set(POS(BQ, lenB - 1, ring), POS(A, lenB - 1, ring), lenB, ring);
        _elem_poly_divrem_basecase(Q, BQ, BQ, 2 * lenB - 1, B, lenB, ring);
        _elem_vec_neg(BQ, BQ, lenB - 1, ring);
        _elem_vec_sub(POS(BQ, lenB - 1, ring), POS(A, lenB - 1, ring), POS(BQ, lenB - 1, ring), lenB, ring);
    }
    else
    {
        const long n2 = lenB / 2;
        const long n1 = lenB - n2;

        elem_ptr W1 = W;
        elem_ptr W2 = POS(W, lenB, ring);

        elem_srcptr p1 = POS(A, 2 * n2, ring);
        elem_srcptr p2;
        elem_srcptr d1 = POS(B, n2, ring);
        elem_srcptr d2 = B;
        elem_srcptr d3 = POS(B, n1, ring);
        elem_srcptr d4 = B;

        elem_ptr q1 = POS(Q, n2, ring);
        elem_ptr q2 = Q;
        elem_ptr dq1 = POS(BQ, n2, ring);
        elem_ptr d1q1 = POS(BQ, 2 * n2, ring);

        elem_ptr d2q1, d3q2, d4q2, t;

        _elem_poly_divrem_divconquer_recursive(q1, d1q1, W1, p1, d1, n1, ring);

        d2q1 = W1;
        _elem_poly_mul(d2q1, q1, n1, d2, n2, ring);

        _elem_vec_swap(dq1, d2q1, n2, ring);
        _elem_vec_add(POS(dq1, n2, ring), POS(dq1, n2, ring), POS(d2q1, n2, ring), n1 - 1, ring);

        t = BQ;
        _elem_vec_sub(t, POS(A, n2 + (n1 - 1), ring), POS(dq1, (n1 - 1), ring), n2, ring);
        p2 = POS(t, -(n2 - 1), ring);

        d3q2 = W1;
        _elem_poly_divrem_divconquer_recursive(q2, d3q2, W2, p2, d3, n2, ring);

        d4q2 = W2;
        _elem_poly_mul(d4q2, d4, n1, q2, n2, ring);

        _elem_vec_swap(BQ, d4q2, n2, ring);
        _elem_vec_add(POS(BQ, n2, ring), POS(BQ, n2, ring), POS(d4q2, n2, ring), n1 - 1, ring);
        _elem_vec_add(POS(BQ, n1, ring), POS(BQ, n1, ring), d3q2, 2 * n2 - 1, ring);
    }
}

static void
__elem_poly_divrem_divconquer(elem_ptr Q, elem_ptr R, elem_srcptr A, long lenA, elem_srcptr B, long lenB, const ring_t ring)
{
    if (lenA < 2 * lenB - 1)
    {
        const long n1 = lenA - lenB + 1;
        const long n2 = lenB - n1;

        elem_srcptr p1 = POS(A, n2, ring);
        elem_srcptr d1 = POS(B, n2, ring);
        elem_srcptr d2 = B;

        elem_ptr W = _elem_vec_init((2 * n1 - 1) + lenB - 1, ring);

        elem_ptr d1q1 = POS(R, n2, ring);
        elem_ptr d2q1 = POS(W, (2 * n1 - 1), ring);

        _elem_poly_divrem_divconquer_recursive(Q, d1q1, W, p1, d1, n1, ring);

        if (n1 >= n2)
            _elem_poly_mul(d2q1, Q, n1, d2, n2, ring);
        else
            _elem_poly_mul(d2q1, d2, n2, Q, n1, ring);

        _elem_vec_swap(R, d2q1, n2, ring);
        _elem_vec_add(POS(R, n2, ring), POS(R, n2, ring), POS(d2q1, n2, ring), n1 - 1, ring);
        _elem_vec_sub(R, A, R, lenA, ring);
        _elem_vec_clear(W, (2 * n1 - 1) + lenB - 1, ring);
    }
    else /* lenA = 2 * lenB - 1 */
    {
        elem_ptr W = _elem_vec_init(lenA, ring);

        _elem_poly_divrem_divconquer_recursive(Q, R, W, A, B, lenB, ring);
        _elem_vec_sub(R, A, R, lenA, ring);

        _elem_vec_clear(W, lenA, ring);
    }
}

void _elem_poly_divrem_divconquer(elem_ptr Q, elem_ptr R, elem_srcptr A, long lenA, elem_srcptr B, long lenB, const ring_t ring)
{
    if (R == A)
    {
        elem_ptr R2 = _elem_vec_init(lenA, ring);
        _elem_poly_divrem_divconquer(Q, R2, A, lenA, B, lenB, ring);
        _elem_vec_swap(R, R2, lenA, ring);
        _elem_vec_clear(R2, lenA, ring);
        return;
    }

    if (lenA <= 2 * lenB - 1)
    {
        __elem_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, ring);
    }
    else /* lenA > 2 * lenB - 1 */
    {
        long shift, n = 2 * lenB - 1;
        elem_ptr QB, W;

        _elem_vec_set(R, A, lenA, ring);
        W = _elem_vec_init(2 * n, ring);
        QB = POS(W, n, ring);

        while (lenA >= n)
        {
            shift = lenA - n;
            _elem_poly_divrem_divconquer_recursive(POS(Q, shift, ring), QB, W, POS(R, shift, ring), B, lenB, ring);
            _elem_vec_sub(POS(R, shift, ring), POS(R, shift, ring), QB, n, ring);
            lenA -= lenB;
        }

        if (lenA >= lenB)
        {
            __elem_poly_divrem_divconquer(Q, W, R, lenA, B, lenB, ring);
            _elem_vec_swap(W, R, lenA, ring);
        }

        _elem_vec_clear(W, 2 * n, ring);
    }
}

void
elem_poly_divrem_divconquer(elem_poly_t Q, elem_poly_t R,
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
        elem_poly_divrem_divconquer(tmp, R, A, B, ring);
        elem_swap(Q, tmp, ring);
        elem_clear(tmp, ring);
        return;
    }

    if (R == A || R == B)
    {
        elem_poly_t tmp;
        elem_init(tmp, ring);
        elem_poly_divrem_divconquer(Q, tmp, A, B, ring);
        elem_swap(R, tmp, ring);
        elem_clear(tmp, ring);
        return;
    }

    elem_poly_fit_length(Q, lenA - lenB + 1, ring);
    elem_poly_fit_length(R, lenA, ring);

    _elem_poly_divrem_divconquer(Q->coeffs, R->coeffs, A->coeffs, lenA, B->coeffs, lenB, ring->parent);

    elem_poly_set_length(Q, lenA - lenB + 1, ring);
    elem_poly_set_length(R, lenA, ring);
    elem_poly_normalise(Q, ring);
    elem_poly_normalise(R, ring);
}

