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

/* check if (p, n) = mullow(poly1, len1, poly2, n, n) where len1 > 0, n >= 0 */
static int
_nmod_poly_mullow_classical_check(mp_srcptr p, mp_srcptr poly1, slong len1,
                            mp_srcptr poly2, slong n, nmod_t mod)
{
    slong i, j, bits, log_len, nlimbs, n1;
    mp_limb_t c;

    len1 = FLINT_MIN(len1, n);

    if (n == 0)
	return 1;

    if (n == 1)
        return p[0] == nmod_mul(poly1[0], poly2[0], mod);

    log_len = FLINT_BIT_COUNT(n);
    bits = FLINT_BITS - (slong) mod.norm;
    bits = 2 * bits + log_len;

    if (bits <= FLINT_BITS)
    {
	for (i = 0; i < n; i++)
        {
            c = 0;
	    
	    for (j = 0; j <= FLINT_MIN(i, len1 - 1); j++)
                c += poly1[j] * poly2[i - j];

	    NMOD_RED(c, c, mod);

	    if (p[i] != c)
               return 0;
        }
    } else
    {
        if (bits <= 2 * FLINT_BITS)
            nlimbs = 2;
        else
            nlimbs = 3;

        for (i = 0; i < n; i++)
        {
            n1 = FLINT_MIN(len1 - 1, i);

            c = _nmod_vec_dot_rev(poly1,
                                  poly2 + i - n1,
                                  n1 + 1, mod, nlimbs);

	    if (p[i] != c)
	       return 0;
        }
    }

    return 1;
}

int _nmod_poly_divides_classical(mp_ptr Q, mp_srcptr A, slong lenA, 
                                           mp_srcptr B, slong lenB, nmod_t mod)
{    
    slong lenQ = lenA - lenB + 1;
    mp_ptr W;
    int res;
    TMP_INIT;

    TMP_START;

    W = TMP_ALLOC(NMOD_DIV_BC_ITCH(lenA, lenB, mod)*sizeof(mp_limb_t));

    _nmod_poly_div_basecase(Q, W, A, lenA, B, lenB, mod);

    /* check coefficients of product one at a time */
    res = _nmod_poly_mullow_classical_check(A, Q, lenQ, B, lenB - 1, mod);
		    
    if (!res)
        _nmod_vec_zero(Q, lenQ);

    TMP_END;

    return res;
}

int nmod_poly_divides_classical(nmod_poly_t Q, const nmod_poly_t A,
		                                           const nmod_poly_t B)
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

    res = _nmod_poly_divides_classical(q, A->coeffs, lenA,
		                                      B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(tQ, Q);
        nmod_poly_clear(tQ);
    }
    
    Q->length = lenA - lenB + 1;
    _nmod_poly_normalise(Q);

    return res;
}
