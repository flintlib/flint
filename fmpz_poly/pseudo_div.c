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

void _fmpz_poly_pseudo_div(fmpz * Q, ulong * d, const fmpz * A, len_t lenA, 
                                                const fmpz * B, len_t lenB)
{
    fmpz * R = _fmpz_vec_init(lenA);
    _fmpz_poly_pseudo_divrem(Q, R, d, A, lenA, B, lenB);
    _fmpz_vec_clear(R, lenA);
}

void fmpz_poly_pseudo_div(fmpz_poly_t Q, ulong * d, const fmpz_poly_t A, 
                                                    const fmpz_poly_t B)
{
    len_t lenq;
    fmpz * q;

    if (B->length == 0)
    {
        printf("Exception (fmpz_poly_pseudo_div). Division by zero.\n");
        abort();
    }
    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        *d = 0;
        return;
    }

    lenq = A->length - B->length + 1;
    if (Q == A || Q == B)
        q = _fmpz_vec_init(lenq);
    else
    {
        fmpz_poly_fit_length(Q, lenq);
        q = Q->coeffs;
    }

    _fmpz_poly_pseudo_div(q, d, A->coeffs, A->length, 
                                B->coeffs, B->length);

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenq;
        Q->length = lenq;
    }
    else
        _fmpz_poly_set_length(Q, lenq);
}
