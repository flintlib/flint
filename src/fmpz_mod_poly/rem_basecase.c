/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_rem_basecase(fmpz *R,
    const fmpz *A, slong lenA, const fmpz *B, slong lenB,
    const fmpz_t invB, const fmpz_mod_ctx_t ctx)
{
    fmpz_t q;
    slong iR;
    fmpz * W;
    TMP_INIT;

    TMP_START;

    fmpz_init(q);

    if (R != A)
    {
        /* cannot use R as it might not have enough space */
        FMPZ_VEC_TMP_INIT(W, lenA);
        _fmpz_vec_set(W, A, lenA);
    }
    else
    {
        W = R;
    }

    for (iR = lenA - 1; iR >= lenB - 1; iR--)
    {
        if (!fmpz_is_zero(W + iR))
        {
            fmpz_mul(q, W + iR, invB);
            fmpz_mod_set_fmpz(q, q, ctx);

            _fmpz_vec_scalar_submul_fmpz(W + (iR - lenB + 1), B, lenB, q);
            _fmpz_mod_vec_set_fmpz_vec(W + (iR - lenB + 1), W + (iR - lenB + 1), lenB, ctx);
        }
    }

    if (R != A)
    {
       _fmpz_vec_set(R, W, lenB - 1);
       FMPZ_VEC_TMP_CLEAR(W, lenA);
    }

    fmpz_clear(q);

    TMP_END;
}

void fmpz_mod_poly_rem_basecase(fmpz_mod_poly_t R, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length;
    fmpz *r;
    fmpz_t invB;

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A, ctx);
        return;
    }

    fmpz_init(invB);
    fmpz_invmod(invB, B->coeffs + (lenB - 1), fmpz_mod_ctx_modulus(ctx));

    if (R == B)
    {
        r = _fmpz_vec_init(lenB - 1);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    _fmpz_mod_poly_rem_basecase(r, A->coeffs, lenA,
                             B->coeffs, lenB, invB, ctx);

    if (R == B)
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
