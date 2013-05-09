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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void _fmpq_poly_lcm(fmpz *L, fmpz_t denL, 
                    const fmpz *A, len_t lenA, const fmpz *B, len_t lenB)
{
    if (lenA == 1)  /* lenA == lenB == 1 */
    {
        fmpz_one(L);
        fmpz_one(denL);
    }
    else
    {
        fmpz *copyA, *copyB;
        fmpz_t s, t;
        len_t lenL;

        fmpz_init(s);
        fmpz_init(t);

        _fmpz_vec_content(s, A, lenA);
        _fmpz_vec_content(t, B, lenB);

        if (fmpz_is_one(s))
        {
            if (fmpz_is_one(t))
            {
                copyA = (fmpz *) A;
                copyB = (fmpz *) B;
            }
            else
            {
                copyA = (fmpz *) A;
                copyB = _fmpz_vec_init(lenB);
                _fmpz_vec_scalar_divexact_fmpz(copyB, B, lenB, t);
            }
        }
        else
        {
            if (fmpz_is_one(s))
            {
                copyA = _fmpz_vec_init(lenA);
                copyB = (fmpz *) B;
                _fmpz_vec_scalar_divexact_fmpz(copyA, A, lenA, s);
            }
            else
            {
                copyA = _fmpz_vec_init(lenA + lenB);
                copyB = copyA + lenA;
                _fmpz_vec_scalar_divexact_fmpz(copyA, A, lenA, s);
                _fmpz_vec_scalar_divexact_fmpz(copyB, B, lenB, t);
            }
        }

        _fmpz_poly_lcm(L, copyA, lenA, copyB, lenB);

        for (lenL = lenA + lenB - 2; !L[lenL]; lenL--) ;
        lenL++;

        fmpz_set(denL, L + (lenL - 1));

        if (A != copyA)
            _fmpz_vec_clear(copyA, lenA + (B != copyB) * lenB);
        else if (B != copyB)
            _fmpz_vec_clear(copyB, lenB);

        fmpz_clear(s);
        fmpz_clear(t);
    }
}

void fmpq_poly_lcm(fmpq_poly_t L, const fmpq_poly_t A, const fmpq_poly_t B)
{
    len_t lenA = A->length, lenB = B->length, lenL = lenA + lenB - 1;

    if (lenA == 0 || lenB == 0)
    {
        fmpq_poly_zero(L);
        return;
    }

    if (L == A || L == B)
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, lenL);
        if (lenA >= lenB)
            _fmpq_poly_lcm(t->coeffs, t->den, A->coeffs, A->length, 
                                              B->coeffs, B->length);
        else
            _fmpq_poly_lcm(t->coeffs, t->den, B->coeffs, B->length, 
                                              A->coeffs, A->length);
        fmpq_poly_swap(t, L);
        fmpq_poly_clear(t);
    }
    else
    {
        fmpq_poly_fit_length(L, lenL);
        if (lenA >= lenB)
            _fmpq_poly_lcm(L->coeffs, L->den, A->coeffs, A->length, 
                                              B->coeffs, B->length);
        else
            _fmpq_poly_lcm(L->coeffs, L->den, B->coeffs, B->length, 
                                              A->coeffs, A->length);
    }

    _fmpq_poly_set_length(L, lenL);
    _fmpq_poly_normalise(L);
}

