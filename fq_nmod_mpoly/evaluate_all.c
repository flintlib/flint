/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void _fq_nmod_mpoly_evaluate_all_fq_nmod_sp(fq_nmod_t ev, const fq_nmod_mpoly_t A,
                  fq_nmod_struct * const * vals, const fq_nmod_mpoly_ctx_t ctx)
{
    ulong l;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong shift, off;
    ulong * ormask, * masks;
    slong * offs;
    slong entries, k_len;
    slong Alen;
    const fq_nmod_struct * Acoeff;
    ulong * Aexp;
    flint_bitcnt_t bits;
    fq_nmod_struct * powers;
    fq_nmod_t t;
    TMP_INIT;

    fq_nmod_init(t, ctx->fqctx);

    Alen = A->length;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    bits = A->bits;

    FLINT_ASSERT(Alen != 0);
    FLINT_ASSERT(bits <= FLINT_BITS);

    TMP_START;

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);

    /* get a mask of present exponent bits */
    ormask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    for (j = 0; j < N; j++)
    {
        ormask[j] = 0;
    }
    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < N; j++)
        {
            ormask[j] |= Aexp[N*i + j];
        }
    }

    /* quick upper bound on number of masks needed */
    entries = N*FLINT_BITS;
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(ulong));
    powers = (fq_nmod_struct *) TMP_ALLOC(entries*sizeof(fq_nmod_struct));

    /* store bit masks for needed powers of two of the variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        FLINT_ASSERT(k < entries);
        mpoly_gen_offset_shift_sp(&off, &shift, i, bits, ctx->minfo);
        fq_nmod_set(t, vals[i], ctx->fqctx);
        for (l = 0; l < bits; l++)
        {
            masks[k] = UWORD(1) << (shift + l);
            if ((masks[k] & ormask[off]) != UWORD(0))
            {
                offs[k] = off;
                fq_nmod_init(powers + k, ctx->fqctx);
                fq_nmod_set(powers + k, t, ctx->fqctx);
                k++;
            }
            fq_nmod_mul(t, t, t, ctx->fqctx);
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len <= entries);

    /* accumulate the final answer */
    fq_nmod_zero(ev, ctx->fqctx);
    for (i = 0; i < Alen; i++)
    {
        fq_nmod_set(t, Acoeff + i, ctx->fqctx);
        for (k = 0; k < k_len; k++)
        {
            if ((Aexp[N*i + offs[k]] & masks[k]) != UWORD(0))
            {
                fq_nmod_mul(t, t, powers + k, ctx->fqctx);
            }
        }
        fq_nmod_add(ev, ev, t, ctx->fqctx);
    }

    for (i = 0; i < k_len; i++)
        fq_nmod_clear(powers + i, ctx->fqctx);

    TMP_END;

    fq_nmod_clear(t, ctx->fqctx);
}

void _fq_nmod_mpoly_evaluate_all_fq_nmod_mp(fq_nmod_t ev, const fq_nmod_mpoly_t A,
                  fq_nmod_struct * const * vals, const fq_nmod_mpoly_ctx_t ctx)
{
    ulong l;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong off;
    ulong * ormask, * masks;
    slong * offs;
    slong entries, k_len;
    slong Alen;
    fq_nmod_struct * Acoeff;
    ulong * Aexp;
    flint_bitcnt_t bits;
    fq_nmod_struct * powers;
    fq_nmod_t t;
    TMP_INIT;

    fq_nmod_init(t, ctx->fqctx);

    Alen = A->length;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    bits = A->bits;

    FLINT_ASSERT(Alen != 0);
    FLINT_ASSERT(bits > FLINT_BITS);

    TMP_START;

    N = mpoly_words_per_exp_mp(bits, ctx->minfo);

    /* get a mask of present exponent bits */
    ormask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    for (j = 0; j < N; j++)
    {
        ormask[j] = 0;
    }
    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < N; j++)
        {
            ormask[j] |= Aexp[N*i + j];
        }
    }

    /* quick upper bound on number of masks needed */
    entries = N*FLINT_BITS;
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(ulong));
    powers = (fq_nmod_struct *) TMP_ALLOC(entries*sizeof(fq_nmod_struct));

    /* store bit masks for needed powers of two of the variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        FLINT_ASSERT(k < entries);
        
        off = mpoly_gen_offset_mp(i, bits, ctx->minfo);
        fq_nmod_set(t, vals[i], ctx->fqctx);
        for (l = 0; l < bits; l++)
        {
            ulong l1 = l/FLINT_BITS;
            ulong l2 = l%FLINT_BITS;
            masks[k] = UWORD(1) << l2;
            if ((masks[k] & ormask[off + l1]) != UWORD(0))
            {
                offs[k] = off + l1;
                fq_nmod_init(powers + k, ctx->fqctx);
                fq_nmod_set(powers + k, t, ctx->fqctx);
                k++;
            }
            fq_nmod_mul(t, t, t, ctx->fqctx);
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len <= entries);

    /* accumulate the final answer */
    fq_nmod_zero(ev, ctx->fqctx);
    for (i = 0; i < Alen; i++)
    {
        fq_nmod_set(t, Acoeff + i, ctx->fqctx);
        for (k = 0; k < k_len; k++)
        {
            if ((Aexp[N*i + offs[k]] & masks[k]) != UWORD(0))
            {
                fq_nmod_mul(t, t, powers + k, ctx->fqctx);
            }
        }
        fq_nmod_add(ev, ev, t, ctx->fqctx);
    }

    for (i = 0; i < k_len; i++)
        fq_nmod_clear(powers + i, ctx->fqctx);

    TMP_END;

    fq_nmod_clear(t, ctx->fqctx);
}


void fq_nmod_mpoly_evaluate_all_fq_nmod(fq_nmod_t ev, const fq_nmod_mpoly_t A,
                  fq_nmod_struct * const * vals, const fq_nmod_mpoly_ctx_t ctx)
{
    if (A->length == 0)
    {
        fq_nmod_zero(ev, ctx->fqctx);
        return;
    }

    if (A->bits <= FLINT_BITS)
    {
        _fq_nmod_mpoly_evaluate_all_fq_nmod_sp(ev, A, vals, ctx);
    }
    else
    {
        _fq_nmod_mpoly_evaluate_all_fq_nmod_mp(ev, A, vals, ctx);
    }
}

