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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

static void
__fmpz_poly_div_divconquer(fmpz * Q, const fmpz * A, long lenA, 
                                     const fmpz * B, long lenB)
{
    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 q1 - 1 by q1 division
         */

        const long q1 = lenA - lenB + 1;
        const long q2 = lenB - q1;

        fmpz * temp = _fmpz_vec_init(2 * q1 - 1);

        _fmpz_poly_div_divconquer_recursive(Q, temp, A + q2, B + q2, q1);

        _fmpz_vec_clear(temp, 2 * q1 - 1);
    }
    else  /* lenA = 2 lenB - 1 */
    {
        fmpz * temp = _fmpz_vec_init(lenA);

        _fmpz_poly_div_divconquer_recursive(Q, temp, A, B, lenB);

        _fmpz_vec_clear(temp, lenA);
    }
}

/* needed due to partial overlap */
static void
_fmpz_vec_sub_dec(fmpz * a, const fmpz * b, const fmpz * c, long n)
{
    long i;

    for (i = n - 1; i >= 0; i--)
        fmpz_sub(a + i, b + i, c + i);
}

void _fmpz_poly_div_divconquer(fmpz *Q, 
                               const fmpz *A, long lenA, 
                               const fmpz *B, long lenB)
{
    if (lenA <= 2 * lenB - 1)
    {
        __fmpz_poly_div_divconquer(Q, A, lenA, B, lenB);
    }
    else  /* lenA > 2 * lenB - 1 */
    {
        fmpz *S, *T;
        long shift, next, n = 2 * lenB - 1;

        S = _fmpz_vec_init(2 * n);
        T = S + n;

        /* To avoid copying all of A, we let S be a window of the
           remainder, taking up to n coefficients at a time */
        shift = lenA - n;
        _fmpz_vec_set(S, A + shift, n);

        while (lenA >= n)
        {
            shift = lenA - n;
            _fmpz_poly_divremlow_divconquer_recursive(Q + shift, T, S, B, lenB);
            next = FLINT_MIN(lenB, shift);
            _fmpz_vec_sub_dec(S + next, S, T, lenB - 1);
            _fmpz_vec_set(S, A + shift - next, next);
            lenA -= lenB;
        }

        if (lenA >= lenB)
            __fmpz_poly_div_divconquer(Q, S, lenA, B, lenB);

        _fmpz_vec_clear(S, 2 * n);
    }
}

void
fmpz_poly_div_divconquer(fmpz_poly_t Q, 
                         const fmpz_poly_t A, const fmpz_poly_t B)
{
    const long lenA = A->length;
    const long lenB = B->length;
    fmpz_poly_t t;
    fmpz *q;

    if (lenB == 0)
    {
        printf("Exception (fmpz_poly_div_divconquer). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_poly_init2(t, lenA - lenB + 1);
        q = t->coeffs;
    }
    else
    {
        fmpz_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    _fmpz_poly_div_divconquer(q, A->coeffs, lenA, B->coeffs, lenB);

    if (Q == A || Q == B)
    {
        _fmpz_poly_set_length(t, lenA - lenB + 1);
        fmpz_poly_swap(t, Q);
        fmpz_poly_clear(t);
    }
    else
        _fmpz_poly_set_length(Q, lenA - lenB + 1);

    _fmpz_poly_normalise(Q);
}
