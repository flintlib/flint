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

    Copyright (C) 2013 William Hart
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_divrem_preinv(fmpz * Q, fmpz * A, slong len1, 
                              const fmpz * B, const fmpz * B_inv, slong len2)
{
   slong n = len1 - len2 + 1;
   fmpz * P = _fmpz_vec_init(len2 - 1);

   _fmpz_poly_div_preinv(Q, A, len1, B, B_inv, len2);
   
   if (len2 - 1 > n)
      _fmpz_poly_mullow(P, B, len2 - 1, Q, n, len2 - 1);
   else
      _fmpz_poly_mullow(P, Q, n, B, len2 - 1, len2 - 1);

   _fmpz_poly_sub(A, A, len2 - 1, P, len2 - 1);

   _fmpz_vec_clear(P, len2 - 1);
}

void
fmpz_poly_divrem_preinv(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, 
                                  const fmpz_poly_t B, const fmpz_poly_t B_inv)
{
    fmpz_poly_t tQ, tR;
    fmpz *q, *r;
    slong len1 = A->length, len2 = B->length;
    slong qlen = len1 - len2 + 1;

    if (len1 < len2)
    {
        fmpz_poly_zero(Q);
        fmpz_poly_set(R, A);
        return;
    }

    if (Q == A || Q == B || Q == B_inv)
    {
        fmpz_poly_init2(tQ, qlen);
        q = tQ->coeffs;
    }
    else
    {
        fmpz_poly_fit_length(Q, qlen);
        q = Q->coeffs;
    }

    if (R == B || R == B_inv)
    {
        fmpz_poly_init2(tR, len1);
        r = tR->coeffs;
    }
    else
    {
        fmpz_poly_fit_length(R, len1);
        r = R->coeffs;
    }

    if (R == B || R == B_inv || R != A)
       _fmpz_vec_set(r, A->coeffs, A->length);

    _fmpz_poly_divrem_preinv(q, r, len1, B->coeffs, B_inv->coeffs, len2);

    if (Q == A || Q == B || Q == B_inv)
    {
        _fmpz_poly_set_length(tQ, qlen);
        fmpz_poly_swap(tQ, Q);
        fmpz_poly_clear(tQ);
    }
    else
        _fmpz_poly_set_length(Q, qlen);

    if (R == B || R == B_inv)
    {
        _fmpz_poly_set_length(tR, len2 - 1);
        fmpz_poly_swap(tR, R);
        fmpz_poly_clear(tR);
    }
    else
        _fmpz_poly_set_length(R, len2 - 1);

    /* no need to normalise Q */
    _fmpz_poly_normalise(R);
}
