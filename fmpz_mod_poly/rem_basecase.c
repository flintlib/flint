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

void _fmpz_mod_poly_rem_basecase(fmpz *R, 
    const fmpz *A, len_t lenA, const fmpz *B, len_t lenB, 
    const fmpz_t invB, const fmpz_t p)
{
    fmpz_t q;
    len_t iR;

    fmpz_init(q);

    if (R != A)
        _fmpz_vec_set(R, A, lenA);

    for (iR = lenA - 1; iR >= lenB - 1; iR--)
    {
        if (!fmpz_is_zero(R + iR))
        {
            fmpz_mul(q, R + iR, invB);
            fmpz_mod(q, q, p);

            _fmpz_vec_scalar_submul_fmpz(R + (iR - lenB + 1), B, lenB, q);
            _fmpz_vec_scalar_mod_fmpz(R + (iR - lenB + 1), R + (iR - lenB + 1), lenB, p);
        }
    }
    fmpz_clear(q);
}

void fmpz_mod_poly_rem_basecase(fmpz_mod_poly_t R, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    const len_t lenA = A->length, lenB = B->length;
    fmpz *r;
    fmpz_t invB;

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A);
        return;
    }

    fmpz_init(invB);
    fmpz_invmod(invB, B->coeffs + (lenB - 1), &(B->p));

    if (R == B)
    {
        r = _fmpz_vec_init(lenA);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenA);
        r = R->coeffs;
    }

    _fmpz_mod_poly_rem_basecase(r, A->coeffs, lenA,
                                   B->coeffs, lenB, invB, &(B->p));

    if (R == B)
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

