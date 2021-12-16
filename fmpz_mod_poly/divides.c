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
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

int _fmpz_mod_poly_divides(fmpz * Q, const fmpz * A, slong lenA, 
                          const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
{    
    fmpz * R;
    fmpz_t invB;
    slong i, lenQ = lenA - lenB + 1;
    int res = 1;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);

    if (lenA < 40 && lenB < 20) 
        return _fmpz_mod_poly_divides_classical(Q, A, lenA, B, lenB, ctx);

    R = _fmpz_vec_init(lenB - 1);

    fmpz_init(invB);
    fmpz_invmod(invB, B + lenB - 1, p);

    if (lenA < 2*lenB - 1)
    {
	slong offset = 0;
        fmpz * P;

	P = (fmpz *) _fmpz_vec_init(2*lenQ - 1);

	_fmpz_vec_zero(R, lenB - 1);

	_fmpz_mod_poly_div(Q, A, lenA, B, lenB, invB, p);

        while (offset < lenB - 1)
        {
            if (offset + 2*lenQ - 1 < lenB)
            {
                _fmpz_mod_poly_mul(P, B + offset, lenQ, Q, lenQ, p);
	        _fmpz_mod_poly_add(R + offset, R + offset, 2*lenQ - 1, P, 2*lenQ - 1, p);
            } else
            {
                _fmpz_mod_poly_mullow(P, Q, lenQ, B + offset, lenQ, p, lenB - offset - 1);
                _fmpz_mod_poly_add(R + offset, R + offset, lenB - offset - 1, P, lenB - offset - 1, p);
            }

	    for (i = 0; i < FLINT_MIN(lenQ, lenB - offset - 1); i++)
            {
                if (!fmpz_equal(R + offset + i, A + offset + i))
		{
                    res = 0;
		    break;
	        }
            }

	    offset += lenQ;
	}

        _fmpz_vec_clear(P, 2*lenQ - 1);
    } else
    {
        _fmpz_mod_poly_divrem(Q, R, A, lenA, B, lenB, invB, p);

        for (i = 0; i < lenB - 1; i++)
        {
            if (!fmpz_is_zero(R + i))
            {
                res = 0;
	        break;
            }
        }
    }

    fmpz_clear(invB);

    _fmpz_vec_clear(R, lenB - 1);

    if (res == 0)
	_fmpz_vec_zero(Q, lenQ);

    return res;
}

int fmpz_mod_poly_divides(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_t tQ;
    fmpz * q;
    slong lenA, lenB;
    int res;

    lenA = A->length;
    lenB = B->length;

    if (lenB == 0 || lenA < lenB)
    {
        fmpz_mod_poly_zero(Q, ctx);
	return lenA == 0;
    }

    if (Q == A || Q == B)
    {
        fmpz_mod_poly_init2(tQ, lenA - lenB + 1, ctx);
        q = tQ->coeffs;
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenA - lenB + 1, ctx);
        q = Q->coeffs;
    }

    res = _fmpz_mod_poly_divides(q, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (Q == A || Q == B)
    {
        fmpz_mod_poly_swap(tQ, Q, ctx);
        fmpz_mod_poly_clear(tQ, ctx);
    }
    
    _fmpz_mod_poly_set_length(Q, lenA - lenB + 1);
    _fmpz_mod_poly_normalise(Q);

    return res;
}
