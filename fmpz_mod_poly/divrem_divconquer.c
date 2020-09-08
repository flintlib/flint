/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

static void
__fmpz_mod_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p)
{
    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 n1 - 1 by n1 division
         */

        const slong n1 = lenA - lenB + 1;
        const slong n2 = lenB - n1;

        const fmpz * p1 = A + n2;
        const fmpz * d1 = B + n2;
        const fmpz * d2 = B;

        fmpz * W = _fmpz_vec_init((2 * n1 - 1) + lenB - 1);

        fmpz * d1q1 = R + n2;
        fmpz * d2q1 = W + (2 * n1 - 1);

        _fmpz_mod_poly_divrem_divconquer_recursive(Q, d1q1, W, p1, d1, n1, 
                                                   invB, p);

        /*
           Compute d2q1 = Q d2, of length lenB - 1
         */

        if (n1 >= n2)
            _fmpz_mod_poly_mul(d2q1, Q, n1, d2, n2, p);
        else
            _fmpz_mod_poly_mul(d2q1, d2, n2, Q, n1, p);

        /*
           Compute BQ = d1q1 * x^n1 + d2q1, of length lenB - 1; 
           then compute R = A - BQ
         */

        _fmpz_vec_swap(R, d2q1, n2);
        _fmpz_mod_poly_add(R + n2, R + n2, n1 - 1, d2q1 + n2, n1 - 1, p);
        _fmpz_mod_poly_sub(R, A, lenA, R, lenA, p);

        _fmpz_vec_clear(W, (2 * n1 - 1) + lenB - 1);
    }
    else  /* lenA = 2 * lenB - 1 */
    {
        fmpz * W = _fmpz_vec_init(lenA);

        _fmpz_mod_poly_divrem_divconquer_recursive(Q, R, W, 
                                                   A, B, lenB, invB, p);

        _fmpz_mod_poly_sub(R, A, lenB - 1, R, lenB - 1, p);

        _fmpz_vec_clear(W, lenA);
    }
}

void _fmpz_mod_poly_divrem_divconquer(fmpz *Q, fmpz *R, 
    const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
    const fmpz_t invB, const fmpz_t p)
{
    if (lenA <= 2 * lenB - 1)
    {
        fmpz * W = _fmpz_vec_init(lenA);
        
        __fmpz_mod_poly_divrem_divconquer(Q, W, A, lenA, B, lenB, invB, p);

        _fmpz_vec_set(R, W, lenB - 1);
        _fmpz_vec_clear(W, lenA);
    }
    else  /* lenA > 2 * lenB - 1 */
    {
        slong shift, n = 2 * lenB - 1, len1;
        fmpz *QB, *W, *S;

        len1 = 2 * n + lenA;
        W = _fmpz_vec_init(len1);
        S = W + 2*n;
        _fmpz_vec_set(S, A, lenA);
        QB = W + n;

        while (lenA >= n)
        {
            shift = lenA - n;
            _fmpz_mod_poly_divrem_divconquer_recursive(Q + shift, QB,
                W, S + shift, B, lenB, invB, p);
            _fmpz_mod_poly_sub(S + shift, S + shift, n, QB, n, p);
            lenA -= lenB;
        }

        if (lenA >= lenB)
        {
            __fmpz_mod_poly_divrem_divconquer(Q, W, S, lenA, B, lenB, invB, p);
            _fmpz_vec_swap(W, S, lenA);
        }

        _fmpz_vec_set(R, S, lenB - 1);
        _fmpz_vec_clear(W, len1);
    }
}

void
fmpz_mod_poly_divrem_divconquer(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                            const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length;
    const slong lenB = B->length;
    const slong lenQ = lenA - lenB + 1;

    fmpz *q, *r;
    fmpz_t invB;

    if (lenB == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A, ctx);
            fmpz_mod_poly_zero(R, ctx);
            return;
        }
        else
        {
            flint_printf("Exception (fmpz_mod_poly_div_basecase). Division by zero.\n");
            flint_abort();
        }
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A, ctx);
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

	if (B->length < 8)
	{
        fmpz_mod_poly_divrem_basecase(Q, R, A, B, ctx);
        return;
    }
	
    fmpz_init(invB);
    fmpz_invmod(invB, fmpz_mod_poly_lead(B, ctx), fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = _fmpz_vec_init(lenB - 1);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    _fmpz_mod_poly_divrem_divconquer(q, r, A->coeffs, lenA, 
                             B->coeffs, lenB, invB, fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fmpz_mod_poly_set_length(Q, lenQ);
    }

    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenB - 1;
        R->length = lenB - 1;
    }

    _fmpz_mod_poly_set_length(R, lenB - 1);
    _fmpz_mod_poly_normalise(R);

    fmpz_clear(invB);
}
