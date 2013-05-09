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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_divrem_f(fmpz_t f, fmpz *Q, fmpz *R, 
                             const fmpz *A, len_t lenA, 
                             const fmpz *B, len_t lenB, const fmpz_t p)
{
    fmpz_t invB;

    fmpz_init(invB);
    fmpz_gcdinv(f, invB, B + lenB - 1, p);

    if (fmpz_is_one(f))
    {
        _fmpz_mod_poly_divrem(Q, R, A, lenA, B, lenB, invB, p);
    }

    fmpz_clear(invB);
}

void fmpz_mod_poly_divrem_f(fmpz_t f, fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                            const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    const len_t lenA = A->length;
    const len_t lenB = B->length;
    const len_t lenQ = lenA - lenB + 1;

    fmpz *q, *r;
    fmpz_t invB;

    fmpz_init(invB);
    fmpz_gcdinv(f, invB, fmpz_poly_lead(B), &(B->p));

    if (!fmpz_is_one(f))
    {
        fmpz_clear(invB);
        return;
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A);
        fmpz_mod_poly_zero(Q);
        fmpz_clear(invB);
        return;
    }

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = _fmpz_vec_init(lenA);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenA);
        r = R->coeffs;
    }

    _fmpz_mod_poly_divrem_divconquer(q, r, A->coeffs, lenA, 
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

    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenA;
        R->length = lenA;
    }
    _fmpz_mod_poly_set_length(R, lenB - 1);
    _fmpz_mod_poly_normalise(R);

    fmpz_clear(invB);
}
