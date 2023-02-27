/*
    Copyright (C) 2021 William Hart

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

int _nmod_poly_divides(mp_ptr Q, mp_srcptr A, slong lenA, 
                                           mp_srcptr B, slong lenB, nmod_t mod)
{    
    mp_ptr R;
    slong i, lenQ = lenA - lenB + 1;
    int res = 1;

    if (lenA < 40 && lenB < 20) 
        return _nmod_poly_divides_classical(Q, A, lenA, B, lenB, mod);

    R = _nmod_vec_init(lenB - 1);

    if (lenA < 2*lenB - 1)
    {
	slong offset = 0;
        mp_ptr P;

	P = (mp_ptr) _nmod_vec_init(2*lenQ - 1);

	_nmod_vec_zero(R, lenB - 1);

	_nmod_poly_div(Q, A, lenA, B, lenB, mod);

        while (offset < lenB - 1)
        {
            if (offset + 2*lenQ - 1 < lenB)
            {
                _nmod_poly_mul(P, B + offset, lenQ, Q, lenQ, mod);
	        _nmod_poly_add(R + offset, R + offset, 2*lenQ - 1, P, 2*lenQ - 1, mod);
            } else
            {
                _nmod_poly_mullow(P, Q, lenQ, B + offset, lenQ, lenB - offset - 1, mod);
                _nmod_poly_add(R + offset, R + offset, lenB - offset - 1, P, lenB - offset - 1, mod);
            }

	    for (i = 0; i < FLINT_MIN(lenQ, lenB - offset - 1); i++)
            {
                if (R[offset + i] != A[offset + i])
		{
                    res = 0;
		    break;
	        }
            }

	    offset += lenQ;
	}

        _nmod_vec_clear(P);
    } else
    {
        _nmod_poly_divrem(Q, R, A, lenA, B, lenB, mod);

        for (i = 0; i < lenB - 1; i++)
        {
            if (R[i] != 0)
            {
                res = 0;
	        break;
            }
        }
    }

    _nmod_vec_clear(R);

    if (res == 0)
	_nmod_vec_zero(Q, lenQ);

    return res;
}

int nmod_poly_divides(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tQ;
    mp_ptr q;
    slong lenA, lenB;
    int res;

    lenA = A->length;
    lenB = B->length;

    if (lenB == 0 || lenA < lenB)
    {
        nmod_poly_zero(Q);
	return lenA == 0;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2(tQ, A->mod.n, lenA - lenB + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    res = _nmod_poly_divides(q, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(tQ, Q);
        nmod_poly_clear(tQ);
    }
    
    Q->length = lenA - lenB + 1;
    _nmod_poly_normalise(Q);

    return res;
}
