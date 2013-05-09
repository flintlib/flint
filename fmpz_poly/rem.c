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

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_rem(fmpz * R, const fmpz * A, len_t lenA, const fmpz * B, len_t lenB)
{
    if (lenA < 15)
    {
        _fmpz_poly_rem_basecase(R, A, lenA, B, lenB);
    }
    else
    {
        fmpz *Q = _fmpz_vec_init(lenA - lenB + 1);

        _fmpz_poly_divrem(Q, R, A, lenA, B, lenB);

        _fmpz_vec_clear(Q, lenA - lenB + 1);
    }
}

void fmpz_poly_rem(fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
    const len_t lenA = A->length, lenB = B->length;
    fmpz *r;

    if (lenB == 0)
    {
        printf("Exception (fmpz_poly_rem). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        fmpz_poly_set(R, A);
        return;
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

    _fmpz_poly_rem(r, A->coeffs, lenA, B->coeffs, lenB);

    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenB - 1;
    }
    _fmpz_poly_set_length(R, lenA);
    _fmpz_poly_normalise(R);
}
