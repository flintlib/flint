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
_fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, ulong * d,
                                  const fmpz * A, len_t lenA, const fmpz * B,
                                  len_t lenB)
{
    const fmpz * leadB = B + (lenB - 1);
    len_t iQ = lenA - lenB, iR = lenA - 1;
    fmpz_t rem;

    fmpz_init(rem);

    *d = 0;
    _fmpz_vec_zero(Q, lenA - lenB + 1);
    if (R != A)
        _fmpz_vec_set(R, A, lenA);

    while (iR >= lenB - 1)
    {
        fmpz_fdiv_qr(Q + iQ, rem, R + iR, leadB);
        if (!fmpz_is_zero(rem))
        {
            _fmpz_vec_scalar_mul_fmpz(Q, Q, lenA - lenB + 1, leadB);
            fmpz_set(Q + iQ, R + iR);
            _fmpz_vec_scalar_mul_fmpz(R, R, lenA, leadB);
            (*d)++;
        }

        if (lenB > 1)
            _fmpz_vec_scalar_submul_fmpz(R + (iR - lenB + 1), B, lenB - 1, Q + iQ);

        fmpz_zero(R + iR);

        iR--;
        iQ--;
    }

    fmpz_clear(rem);
}

void
fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R,
                                 ulong * d, const fmpz_poly_t A,
                                 const fmpz_poly_t B)
{
    len_t lenq, lenr;
    fmpz *q, *r;

    if (B->length == 0)
    {
        printf("Exception (fmpz_poly_pseudo_divrem_basecase). Division by zero.\n");
        abort();
    }
    if (Q == R)
    {
        printf("Exception (fmpz_poly_pseudo_divrem_basecase). \n"
               "Output arguments Q and R may not be aliased.\n");
        abort();
    }
    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        fmpz_poly_set(R, A);
        *d = 0;
        return;
    }

    lenq = A->length - B->length + 1;
    lenr = A->length;
    if (Q == A || Q == B)
        q = _fmpz_vec_init(lenq);
    else
    {
        fmpz_poly_fit_length(Q, lenq);
        q = Q->coeffs;
    }
    if (R == B)
        r = _fmpz_vec_init(lenr);
    else
    {
        fmpz_poly_fit_length(R, lenr);
        r = R->coeffs;
    }
    
    _fmpz_poly_pseudo_divrem_basecase(q, r, d, A->coeffs, A->length, 
                                               B->coeffs, B->length);
    
    for (lenr = B->length - 2; (lenr >= 0) && !r[lenr]; lenr--) ;
    lenr++;
    
    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc = lenq;
        Q->length = lenq;
    }
    else
        _fmpz_poly_set_length(Q, lenq);
    if (R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc = A->length;
        R->length = lenr;
    }
    else
        _fmpz_poly_set_length(R, lenr);
}
