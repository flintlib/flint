/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
    leading coefficient wrt gen(0), ..., gen(num_vars-1)
    c will be returned with c->bits == A->bits
*/
void fq_nmod_mpolyl_lead_coeff(
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_t A,
    slong num_vars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d;
    slong i, j, off, shift;
    ulong mask, first_mask;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexps = A->exps;
    ulong * cexps;
    slong Alen = A->length;

    FLINT_ASSERT(c != A);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(0 < num_vars && num_vars < ctx->minfo->nvars);

    mpoly_gen_offset_shift_sp(&off, &shift, num_vars - 1, A->bits, ctx->minfo);

    mask = (-UWORD(1)) << shift;

    i = 0;
    first_mask = (Aexps + N*i)[off] & mask;

    for (i = 1; i < Alen; i++)
    {
        if (((Aexps + N*i)[off] & mask) != first_mask)
            goto break_outer;

        for (j = off + 1; j < N; j++)
            if ((Aexps + N*(i - 1))[j] != (Aexps + N*i)[j])
                goto break_outer;
    }

break_outer:

    fq_nmod_mpoly_fit_length_reset_bits(c, i, A->bits, ctx);
    c->length = i;
    cexps = c->exps;

    d = fq_nmod_ctx_degree(ctx->fqctx);
    _nmod_vec_set(c->coeffs, A->coeffs, d*i);

    mask = ~mask;

    for (i = 0; i < c->length; i++)
    {
        for (j = 0; j < off; j++)
            (cexps + N*i)[j] = (Aexps + N*i)[j];

        (cexps + N*i)[off] = mask & (Aexps + N*i)[off];

        for (j = off + 1; j < N; j++)
            (cexps + N*i)[j] = 0;
    }
}

