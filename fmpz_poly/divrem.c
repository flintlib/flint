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

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_divrem(fmpz * Q, fmpz * R,
                       const fmpz * A, long lenA, const fmpz * B, long lenB)
{
    if (lenB < 6)
        _fmpz_poly_divrem_basecase(Q, R, A, lenA, B, lenB);
    else
        _fmpz_poly_divrem_divconquer(Q, R, A, lenA, B, lenB);
}

void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R,
                      const fmpz_poly_t A, const fmpz_poly_t B)
{
    const long lenA = A->length, lenB = B->length;
    fmpz *q, *r;

    if (lenB == 0)
    {
        printf("Exception (fmpz_poly_divrem). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        fmpz_poly_set(R, A);
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenA - lenB + 1);
    }
    else
    {
        fmpz_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = _fmpz_vec_init(lenA);
    }
    else
    {
        fmpz_poly_fit_length(R, lenA);
        r = R->coeffs;
    }

    _fmpz_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB);

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenA - lenB + 1;
    }
    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenB - 1;
    }
    _fmpz_poly_set_length(Q, lenA - lenB + 1);
    _fmpz_poly_set_length(R, lenA);
    _fmpz_poly_normalise(Q);
    _fmpz_poly_normalise(R);
}
