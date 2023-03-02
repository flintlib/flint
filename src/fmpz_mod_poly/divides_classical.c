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

/* check if (p, n) = mullow(poly1, len1, poly2, n, n) where len1 > 0, n >= 0 */
static int
_fmpz_mod_poly_mullow_classical_check(const fmpz * p,
                            const fmpz * poly1, slong len1,
                        const fmpz * poly2, slong n, const fmpz_mod_ctx_t mod)
{
    slong i, n1;
    fmpz_t c;
    int res = 1;

    len1 = FLINT_MIN(len1, n);

    if (n == 0)
	return 1;

    fmpz_init(c);

    if (n == 1)
    {
        fmpz_mod_mul(c, poly1 + 0, poly2 + 0, mod);

        res = fmpz_equal(p + 0, c);
    } else
    {
        for (i = 0; i < n; i++)
        {
            n1 = FLINT_MIN(len1 - 1, i);

            _fmpz_mod_vec_dot_rev(c, poly1,
                                     poly2 + i - n1,
                                     n1 + 1, mod);

            if (!fmpz_equal(p + i, c))
	    {
                res = 0;
                break;
            }
        }
    }

    fmpz_clear(c);

    return res;
}

int _fmpz_mod_poly_divides_classical(fmpz * Q, const fmpz * A, slong lenA, 
                          const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
{    
    slong lenQ = lenA - lenB + 1;
    fmpz * W;
    int res;
    fmpz_t invB;

    W = _fmpz_vec_init(lenA);

    fmpz_init(invB);
    fmpz_invmod(invB, B + lenB - 1, fmpz_mod_ctx_modulus(ctx));

    _fmpz_mod_poly_div_basecase(Q, W, A, lenA,
		                     B, lenB, invB, fmpz_mod_ctx_modulus(ctx));
    /* check coefficients of product one at a time */
    res = _fmpz_mod_poly_mullow_classical_check(A, Q, lenQ, B, lenB - 1, ctx);
		    
    if (!res)
        _fmpz_vec_zero(Q, lenQ);

    fmpz_clear(invB);
    _fmpz_vec_clear(W, lenA);

    return res;
}

int fmpz_mod_poly_divides_classical(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
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

    res = _fmpz_mod_poly_divides_classical(q, A->coeffs, lenA,
		                                      B->coeffs, lenB, ctx);

    if (Q == A || Q == B)
    {
        fmpz_mod_poly_swap(tQ, Q, ctx);
        fmpz_mod_poly_clear(tQ, ctx);
    }
    
    _fmpz_mod_poly_set_length(Q, lenA - lenB + 1);
    _fmpz_mod_poly_normalise(Q);

    return res;
}
