/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#  if HAVE_ALLOCA_H
#   include <alloca.h>
#  else
#   if _MSC_VER
#    include <malloc.h>
#    define alloca _alloca
#   else
#    ifdef __DECC
#     define alloca(x) __ALLOCA(x)
#    else
#     ifdef BSD
#      include <stdlib.h>
#     else
#      error Could not find alloca
#     endif
#    endif
#   endif
#  endif
# endif
#endif

#include "flint-impl.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_div_basecase(fmpz *Q, fmpz *R, 
    const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
    const fmpz_t invB, const fmpz_t p)
{
    const slong alloc = (R == NULL) ? lenA : 0;
    slong lenR = lenB - 1, iQ;
    TMP_INIT;

    TMP_START;

    if (alloc)
        FMPZ_VEC_TMP_INIT(R, alloc);

    if (R != A)
        _fmpz_vec_set(R + lenR, A + lenR, lenA - lenR);

    for (iQ = lenA - lenB; iQ >= 0; iQ--)
    {
        if (fmpz_is_zero(R + lenA - 1))
        {
            fmpz_zero(Q + iQ);
        }
        else
        {
            fmpz_mul(Q + iQ, R + lenA - 1, invB);
            fmpz_mod(Q + iQ, Q + iQ, p);

            _fmpz_vec_scalar_submul_fmpz(R + lenA - lenR - 1, B, lenR, Q + iQ);
        }
        if (iQ > 0)
            fmpz_mod(R + lenA - 2, R + lenA - 2, p);

        if (lenR - 1 >= iQ)
        {
            B++;
            lenR--;
        }

        lenA--;
    }

    if (alloc)
        FMPZ_VEC_TMP_CLEAR(R, alloc);

    TMP_END;
}

void fmpz_mod_poly_div_basecase(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                            const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    fmpz *q;
    fmpz_t invB;

    if (lenB == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A, ctx);
            return;
        }
        else
            flint_throw(FLINT_DIVZERO, "fmpz_mod_poly_div_basecase\n");
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    fmpz_init(invB);
    fmpz_invmod(invB, B->coeffs + (lenB - 1), fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    _fmpz_mod_poly_div_basecase(q, NULL, A->coeffs, lenA,
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

    fmpz_clear(invB);
}
