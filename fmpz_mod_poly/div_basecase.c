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

    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

******************************************************************************/

#include <stdlib.h>
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_div_basecase(fmpz *Q, fmpz *R, 
    const fmpz *A, len_t lenA, const fmpz *B, len_t lenB, 
    const fmpz_t invB, const fmpz_t p)
{
    const len_t alloc = (R == NULL) ? lenA : 0;
    len_t lenR = lenB - 1, iQ;

    if (alloc)
        R = _fmpz_vec_init(alloc);
    if (R != A)
        _fmpz_vec_set(R + lenR, A + lenR, lenA - lenR);

    for (iQ = lenA - lenB; iQ >= 0; iQ--)
    {
        if (fmpz_is_zero(R + lenA - 1))
        {
            fmpz_zero(Q + iQ);
        }
        else
        {
            fmpz_mul(Q + iQ, R + lenA - 1, invB);
            fmpz_mod(Q + iQ, Q + iQ, p);

            _fmpz_vec_scalar_submul_fmpz(R + lenA - lenR - 1, B, lenR, Q + iQ);
            _fmpz_vec_scalar_mod_fmpz(R + lenA - lenR - 1, R + lenA - lenR - 1, lenR, p);
        }

        if (lenR - 1 >= iQ)
        {
            B++;
            lenR--;
        }

        lenA--;
    }

    if (alloc)
        _fmpz_vec_clear(R, alloc);
}

void fmpz_mod_poly_div_basecase(fmpz_mod_poly_t Q, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    const len_t lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    fmpz *q;
    fmpz_t invB;

    if (lenA < lenB)
    {
        fmpz_mod_poly_zero(Q);
        return;
    }

    fmpz_init(invB);
    fmpz_invmod(invB, B->coeffs + (lenB - 1), &(B->p));

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }

    _fmpz_mod_poly_div_basecase(q, NULL, A->coeffs, lenA,
                                         B->coeffs, lenB, invB, &(B->p));

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fmpz_mod_poly_set_length(Q, lenQ);
    }

    fmpz_clear(invB);
}

