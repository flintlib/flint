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
    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_divrem(mp_ptr Q, mp_ptr R, mp_srcptr A, long lenA, 
                                  mp_srcptr B, long lenB, nmod_t mod)
{
    if (lenB < 6000)
        _nmod_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, mod);
    else
        _nmod_poly_divrem_newton(Q, R, A, lenA, B, lenB, mod);
}

void
nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R,
                 const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tQ, tR;
    mp_ptr q, r;
    long A_len, B_len;

    B_len = B->length;
    
    if (B_len == 0)
    {
        printf("Exception: division by zero in nmod_poly_divrem\n");
        abort();
    }

    A_len = A->length;
    
    if (A_len < B_len)
    {
        nmod_poly_zero(Q);
        nmod_poly_set(R, A);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2(tQ, A->mod.n, A_len - B_len + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, A_len - B_len + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_poly_init2(tR, A->mod.n, B_len - 1);
        r = tR->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, B_len - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem(q, r, A->coeffs, A_len,
                            B->coeffs, B_len, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(tQ, Q);
        nmod_poly_clear(tQ);
    }
    
    Q->length = A_len - B_len + 1;

    if (R == A || R == B)
    {
        nmod_poly_swap(tR, R);
        nmod_poly_clear(tR);
    }
        
    R->length = B_len - 1;

    _nmod_poly_normalise(R);
}
