/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_repack_bits(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                              flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fq_nmod_mpoly_t T;

    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    if (B->bits == Abits || B->length == 0)
    {
        fq_nmod_mpoly_set(A, B, ctx);
        return 1;
    }
    
    /* must use B->alloc because we are going to swap coeff in aliasing case */
    fq_nmod_mpoly_init3(T, B->alloc, Abits, ctx);
    success = mpoly_repack_monomials(T->exps, Abits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    if (success)
    {
        if (A == B)
        {
            fq_nmod_struct * temp = A->coeffs;
            A->coeffs = T->coeffs;
            T->coeffs = temp;
        }
        else
        {
            for (i = 0; i < B->length; i++)
            {
                fq_nmod_set(T->coeffs + i, B->coeffs + i, ctx->fqctx);
            }
        }
        _fq_nmod_mpoly_set_length(T, B->length, ctx);
        fq_nmod_mpoly_swap(A, T, ctx);
    }

    fq_nmod_mpoly_clear(T, ctx);

    return success;
}
