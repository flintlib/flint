/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int nmod_mpoly_repack_bits(nmod_mpoly_t A, const nmod_mpoly_t B,
                                 mp_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx)
{
    int success;
    nmod_mpoly_t T;

    if (B->bits == Abits || B->length == 0)
    {
        nmod_mpoly_set(A, B, ctx);
        return 1;
    }
    
    /* must use B->alloc because we are going to swap coeff in aliasing case */
    nmod_mpoly_init3(T, B->alloc, Abits, ctx);
    success = mpoly_repack_monomials(T->exps, Abits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    if (success)
    {
        if (A == B)
        {
            mp_limb_t * temp = A->coeffs;
            A->coeffs = T->coeffs;
            T->coeffs = temp;
        }
        else
        {
            _nmod_vec_set(T->coeffs, B->coeffs, B->length);
        }
        _nmod_mpoly_set_length(T, B->length, ctx);
        nmod_mpoly_swap(A, T, ctx);
    }

    nmod_mpoly_clear(T, ctx);

    return success;
}
