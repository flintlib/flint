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

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_div_basecase(fmpz * Q, fmpz * R, const fmpz * A, long lenA,
                        const fmpz * B, long lenB)
{
    const fmpz * leadB = B + lenB - 1;
    long B1, iQ = lenA - lenB;
    long alloc;

    while (lenA >= lenB && fmpz_cmpabs(A + lenA - 1, leadB) < 0)
    {
        fmpz_zero(Q + iQ);
        iQ--;
        lenA--;
    }

    if (lenA < lenB)
        return;

    alloc = (R == NULL) ? lenA : 0;
    if (alloc)
        R = _fmpz_vec_init(alloc);
    if (R != A)
        _fmpz_vec_set(R + lenB - 1, A + lenB - 1, lenA - lenB + 1);

    B1 = lenB - 1;
    
    while (lenA >= lenB)
    {
        if (fmpz_cmpabs(R + lenA - 1, leadB) < 0)
            fmpz_zero(Q + iQ);
        else
        {
            fmpz_fdiv_q(Q + iQ, R + lenA - 1, leadB);
            _fmpz_vec_scalar_submul_fmpz(R + lenA - B1 - 1, B, B1, Q + iQ);
        }

        if (B1 >= lenA - lenB + 1)
        {
            B++;
            B1--;
        }

        lenA--;
        iQ--;
    }

    if (alloc)
        _fmpz_vec_clear(R, alloc);
}

void
fmpz_poly_div_basecase(fmpz_poly_t Q, 
                       const fmpz_poly_t A, const fmpz_poly_t B)
{
    long lenq;
    fmpz *q;
    
    if (B->length == 0)
    {
        printf("Exception (fmpz_poly_div_basecase). Division by zero.\n");
        abort();
    }
    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        return;
    }

    lenq = A->length - B->length + 1;
    if ((Q == A) || (Q == B))
        q = _fmpz_vec_init(lenq);
    else
    {
        fmpz_poly_fit_length(Q, lenq);
        q = Q->coeffs;
    }

    _fmpz_poly_div_basecase(q, NULL, A->coeffs, A->length,
                                     B->coeffs, B->length);

    if ((Q == A) || (Q == B))
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc = lenq;
        Q->length = lenq;
    }
    else
        _fmpz_poly_set_length(Q, lenq);

    _fmpz_poly_normalise(Q);
}
