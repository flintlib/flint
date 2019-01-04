/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


int fmpz_mpoly_repack_bits(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                 mp_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpoly_t T;

    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    if (B->bits == Abits || B->length == 0)
    {
        fmpz_mpoly_set(A, B, ctx);
        return 1;
    }
    
    /* must use B->alloc because we are going to swap coeff in aliasing case */
    fmpz_mpoly_init3(T, B->alloc, Abits, ctx);
    success = mpoly_repack_monomials(T->exps, Abits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    if (success)
    {
        if (A == B)
        {
            fmpz * temp = A->coeffs;
            A->coeffs = T->coeffs;
            T->coeffs = temp;
        }
        else
        {
            _fmpz_vec_set(T->coeffs, B->coeffs, B->length);
        }
        _fmpz_mpoly_set_length(T, B->length, ctx);
        fmpz_mpoly_swap(A, T, ctx);
    }

    fmpz_mpoly_clear(T, ctx);

    return success;
}
