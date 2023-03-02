/*
    Copyright (C) 2008, 2009, 2011 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

static void
__fmpz_mod_poly_div_divconquer(fmpz * Q, const fmpz * A, slong lenA, 
                  const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_t p)
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
        
        fmpz * W = (fmpz *) _fmpz_vec_init(2*lenB);

        _fmpz_mod_poly_div_divconquer_recursive(Q, W, p1, d1, n1, invB, p);

        _fmpz_vec_clear(W, 2*lenB);
    }
    else  /* lenA = 2 * lenB - 1 */
    {
        fmpz * W = _fmpz_vec_init(2*lenB);
 
        _fmpz_mod_poly_div_divconquer_recursive(Q, W, A, B, lenB, invB, p);
        
        _fmpz_vec_clear(W, 2*lenB);
    }
}


/* needed due to partial overlap */
static void
_fmpz_mod_poly_sub_dec(fmpz * a, const fmpz * b,
                                        const fmpz * c, slong n, const fmpz * p)
{
    slong i;

    for (i = n - 1; i >= 0; i--)
    {
        fmpz_sub(a + i, b + i, c + i);

        if (fmpz_sgn(a + i) < 0)
            fmpz_add(a + i, a + i, p);
    }
}

void
_fmpz_mod_poly_div_divconquer(fmpz * Q, const fmpz * A, slong lenA,
                  const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_t p)
{
    if (lenA <= 2 * lenB - 1)
    {
        __fmpz_mod_poly_div_divconquer(Q, A, lenA, B, lenB, invB, p);
    }
    else  /* lenA > 2 * lenB - 1 */
    {
        fmpz * S, * T, * R;
        slong shift, next, n = 2 * lenB - 1;

        S = _fmpz_vec_init(3 * n);
        T = S + n;
        R = T + n;

        shift = lenA - n;
        _fmpz_vec_set(S, A + shift, n);

        while (lenA >= n)
        {
            shift = lenA - n;
            _fmpz_mod_poly_divrem_divconquer_recursive(Q + shift, T, R, S, B, lenB, invB, p);
            next = FLINT_MIN(lenB, shift);
            _fmpz_mod_poly_sub_dec(S + next, S, T, lenB - 1, p);
            _fmpz_vec_set(S, A + shift - next, next);
            lenA -= lenB;
        }

        if (lenA >= lenB)
            __fmpz_mod_poly_div_divconquer(Q, S, lenA, B, lenB, invB, p);

        _fmpz_vec_clear(S, 3*n);
    }
}


void
fmpz_mod_poly_div_divconquer(fmpz_mod_poly_t Q,
      const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz * q;
    const slong Alen = A->length;
    const slong Blen = B->length;
    const slong Qlen = Alen - Blen + 1;
    fmpz_t invB;

    if (Blen == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A, ctx);
            return;
        } else
        {
            flint_printf("Exception (fmpz_mod_poly_div_divconquer). Division by zero.\n");
            flint_abort();
        }
    }

    if (Alen < Blen)
    {
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    if (Q->length < 8)
    {
        fmpz_mod_poly_div_basecase(Q, A, B, ctx);
        return;
    }

    fmpz_init(invB);
    fmpz_invmod(invB, fmpz_mod_poly_lead(B, ctx), fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(Qlen);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, Qlen, ctx);
        q = Q->coeffs;
    }

    _fmpz_mod_poly_div_divconquer(q, A->coeffs, Alen,
                              B->coeffs, Blen, invB, fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc = Qlen;
        Q->length = Qlen;
    } else
        _fmpz_mod_poly_set_length(Q, Qlen);
    
    fmpz_clear(invB);
}
