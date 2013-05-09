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

void _fmpz_poly_pseudo_rem(fmpz * R, ulong * d, const fmpz * A, len_t lenA, 
                                                const fmpz * B, len_t lenB)
{
    fmpz * Q = _fmpz_vec_init(lenA + lenB - 1);
    _fmpz_poly_pseudo_divrem(Q, R, d, A, lenA, B, lenB);
    _fmpz_vec_clear(Q, lenA + lenB - 1);
}

void fmpz_poly_pseudo_rem(fmpz_poly_t R, ulong * d, const fmpz_poly_t A, 
                                                    const fmpz_poly_t B)
{
    len_t lenr;
    fmpz * r;

    if (B->length == 0)
    {
        printf("Exception (fmpz_poly_pseudo_rem). Division by zero.\n");
        abort();
    }
    if (A->length < B->length)
    {
        fmpz_poly_set(R, A);
        *d = 0;
        return;
    }

    lenr = A->length;
    if (R == B || R == A)
        r = _fmpz_vec_init(lenr);
    else
    {
        fmpz_poly_fit_length(R, lenr);
        r = R->coeffs;
    }

    _fmpz_poly_pseudo_rem(r, d, A->coeffs, A->length, 
                                B->coeffs, B->length);

    for (lenr = B->length - 2; (lenr >= 0) && !r[lenr]; lenr--) ;
    lenr++;

    if (R == B || R == A)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = A->length;
        R->length = lenr;
    }
    else
        _fmpz_poly_set_length(R, lenr);
}
