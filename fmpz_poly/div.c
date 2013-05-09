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

void
_fmpz_poly_div(fmpz * Q, const fmpz * A, len_t lenA, const fmpz * B, len_t lenB)
{
    _fmpz_poly_div_divconquer(Q, A, lenA, B, lenB);
}

void
fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
    fmpz_poly_t tQ;
    fmpz *q;

    if (B->length == 0)
    {
        printf("Exception (fmpz_poly_div). Division by zero.\n");
        abort();
    }

    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_poly_init2(tQ, A->length - B->length + 1);
        q = tQ->coeffs;
    }
    else
    {
        fmpz_poly_fit_length(Q, A->length - B->length + 1);
        q = Q->coeffs;
    }

    _fmpz_poly_div(q, A->coeffs, A->length, B->coeffs, B->length);

    if (Q == A || Q == B)
    {
        _fmpz_poly_set_length(tQ, A->length - B->length + 1);
        fmpz_poly_swap(tQ, Q);
        fmpz_poly_clear(tQ);
    }
    else
        _fmpz_poly_set_length(Q, A->length - B->length + 1);

    _fmpz_poly_normalise(Q);
}
