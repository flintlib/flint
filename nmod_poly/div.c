/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_div(mp_ptr Q, mp_srcptr A, slong lenA, 
                                  mp_srcptr B, slong lenB, nmod_t mod)
{
    TMP_INIT;
    
    if (lenB < 15)
    {
        mp_ptr W;
        
        TMP_START;
        W = TMP_ALLOC(NMOD_DIV_BC_ITCH(lenA, lenB, mod)*sizeof(mp_limb_t));
        _nmod_poly_div_basecase(Q, W, A, lenA, B, lenB, mod);
        TMP_END;
    }
    else if (lenB < 6000)
        _nmod_poly_div_divconquer(Q, A, lenA, B, lenB, mod);
    else
        _nmod_poly_div_newton(Q, A, lenA, B, lenB, mod);
}

void
nmod_poly_div(nmod_poly_t Q, 
                 const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tQ;
    mp_ptr q;
    slong A_len, B_len;

    B_len = B->length;
    
    if (B_len == 0)
    {
        if (nmod_poly_modulus(B) == 1)
        {
            nmod_poly_set(Q, A);
            return;
        } else
        {                                                                                
            flint_printf("Exception (nmod_poly_divrem). Division by zero.\n");
            flint_abort();
        }
    }

    A_len = A->length;
    
    if (A_len < B_len)
    {
        nmod_poly_zero(Q);
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

    _nmod_poly_div(q, A->coeffs, A_len,
                            B->coeffs, B_len, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(tQ, Q);
        nmod_poly_clear(tQ);
    }
    
    Q->length = A_len - B_len + 1;
}
