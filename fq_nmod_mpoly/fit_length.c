/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void _fq_nmod_mpoly_fit_length(fq_nmod_struct ** coeff,
                              ulong ** exps, slong * alloc, slong len, slong N,
                                                     const fq_nmod_ctx_t fqctx)
{
    if (len > *alloc)
    {
        slong i;
        len = FLINT_MAX(len, 2*(*alloc));
        (* coeff) = (fq_nmod_struct *) flint_realloc(* coeff,
                                                   len*sizeof(fq_nmod_struct));
        (* exps) = (ulong *) flint_realloc(*exps, len*N*sizeof(ulong)); 
        for (i = *alloc; i < len; i++)
            fq_nmod_init((* coeff) + i, fqctx);
        (* alloc) = len;
    }
}

void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, slong length,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
            A->coeffs = (fq_nmod_struct *) flint_malloc(new_alloc
                                                      *sizeof(fq_nmod_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
            A->coeffs = (fq_nmod_struct *) flint_realloc(A->coeffs,
                                             new_alloc*sizeof(fq_nmod_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fq_nmod_init(A->coeffs + i, ctx->fqctx);
        }
        A->alloc = new_alloc;
    }
}
