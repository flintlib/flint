/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"
#include "fq_zech_mpoly_factor.h"

/*
    leading coefficient wrt gen(0), ..., gen(num_vars-1)
    c will be returned with c->bits == A->bits
*/
void fq_zech_mpolyl_lead_coeff(
    fq_zech_mpoly_t c,
    const fq_zech_mpoly_t A,
    slong num_vars,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, j, off, shift;
    ulong old_shift, new_shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexps = A->exps;
    ulong * cexps;
    slong Alen = A->length;

    FLINT_ASSERT(c != A);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(0 < num_vars && num_vars < ctx->minfo->nvars);

    mpoly_gen_offset_shift_sp(&off, &shift, num_vars - 1, A->bits, ctx->minfo);

    i = 0;
    old_shift = (Aexps + N*i)[off] >> shift;

    for (i = 1; i < Alen; old_shift = new_shift, i++)
    {
        new_shift = (Aexps + N*i)[off] >> shift;

        if (new_shift != old_shift)
            break;

        for (j = off + 1; j < N; j++)
            if ((Aexps + N*(i - 1))[j] != (Aexps + N*i)[j])
                break;
    }

    fq_zech_mpoly_fit_length_reset_bits(c, i, A->bits, ctx);
    c->length = i;
    cexps = c->exps;

    for (i = 0; i < c->length; i++)
        fq_zech_set(c->coeffs + i, A->coeffs + i, ctx->fqctx);

    for (i = 0; i < c->length; i++)
    {
        for (j = 0; j < off; j++)
            (cexps + N*i)[j] = (Aexps + N*i)[j];

        (cexps + N*i)[off] = ((-UWORD(1)) >> (FLINT_BITS - shift)) & 
                                                            (Aexps + N*i)[off];
        for (j = off + 1; j < N; j++)
            (cexps + N*i)[j] = 0;
    }
}

