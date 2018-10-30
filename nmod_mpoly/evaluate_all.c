/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


ulong _nmod_mpoly_evaluate_all_ui_sp(nmod_mpoly_t A,
                                const ulong * vals, const nmod_mpoly_ctx_t ctx)
{
    ulong l;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong shift, off;
    ulong * ormask, * masks;
    slong * offs;
    slong entries, k_len;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    mp_bitcnt_t bits;
    mp_limb_t * powers;
    mp_limb_t t, r, acc0, acc1, acc2, pp0, pp1;
    TMP_INIT;

    Alen = A->length;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    bits = A->bits;

    FLINT_ASSERT(Alen != 0);
    FLINT_ASSERT(bits <= FLINT_BITS);

    TMP_START;

    N = mpoly_words_per_exp(bits, ctx->minfo);

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
    powers = (mp_limb_t *) TMP_ALLOC(entries*sizeof(mp_limb_t));

    /* store bit masks for needed powers of two of the variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        FLINT_ASSERT(k < entries);
        mpoly_gen_offset_shift(&off, &shift, i, N, bits, ctx->minfo);
        NMOD_RED(t, vals[i], ctx->ffinfo->mod);
        for (l = 0; l < bits; l++)
        {
            masks[k] = UWORD(1) << (shift + l);
            if ((masks[k] & ormask[off]) != UWORD(0))
            {
                offs[k] = off;
                powers[k] = t;
                k++;
            }
            t = nmod_mul(t, t, ctx->ffinfo->mod);
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len <= entries);

    /* accumulate the final answer */
    acc0 = acc1 = acc2 = 0;    
    for (i = 0; i < Alen; i++)
    {
        t = UWORD(1);
        for (k = 0; k < k_len; k++)
        {
            if ((Aexp[N*i + offs[k]] & masks[k]) != UWORD(0))
            {
                t = nmod_mul(t, powers[k], ctx->ffinfo->mod);
            }
        }
        umul_ppmm(pp1, pp0, Acoeff[i], t);
        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
    }
    NMOD_RED3(r, acc2, acc1, acc0, ctx->ffinfo->mod);

    TMP_END;
    return r;
}



ulong _nmod_mpoly_evaluate_all_ui_mp(nmod_mpoly_t A,
                                const ulong * vals, const nmod_mpoly_ctx_t ctx)
{
    ulong l;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong off;
    ulong * ormask, * masks;
    slong * offs;
    slong entries, k_len;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    mp_bitcnt_t bits;
    mp_limb_t * powers;
    mp_limb_t t, r, acc0, acc1, acc2, pp0, pp1;
    TMP_INIT;

    Alen = A->length;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    bits = A->bits;

    FLINT_ASSERT(Alen != 0);
    FLINT_ASSERT(bits > FLINT_BITS);

    TMP_START;

    N = mpoly_words_per_exp(bits, ctx->minfo);

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
    powers = (mp_limb_t *) TMP_ALLOC(entries*sizeof(mp_limb_t));

    /* store bit masks for needed powers of two of the variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        FLINT_ASSERT(k < entries);
        
        off = mpoly_gen_offset_mp(i, N, bits, ctx->minfo);
        NMOD_RED(t, vals[i], ctx->ffinfo->mod);
        for (l = 0; l < bits; l++)
        {
            ulong l1 = l/FLINT_BITS;
            ulong l2 = l%FLINT_BITS;
            masks[k] = UWORD(1) << l2;
            if ((masks[k] & ormask[off + l1]) != UWORD(0))
            {
                offs[k] = off + l1;
                powers[k] = t;
                k++;
            }
            t = nmod_mul(t, t, ctx->ffinfo->mod);
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len <= entries);

    /* accumulate the final answer */
    acc0 = acc1 = acc2 = 0;    
    for (i = 0; i < Alen; i++)
    {
        t = UWORD(1);
        for (k = 0; k < k_len; k++)
        {
            if ((Aexp[N*i + offs[k]] & masks[k]) != UWORD(0))
            {
                t = nmod_mul(t, powers[k], ctx->ffinfo->mod);
            }
        }
        umul_ppmm(pp1, pp0, Acoeff[i], t);
        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
    }
    NMOD_RED3(r, acc2, acc1, acc0, ctx->ffinfo->mod);

    TMP_END;
    return r;
}


ulong nmod_mpoly_evaluate_all_ui(nmod_mpoly_t A,
                                const ulong * vals, const nmod_mpoly_ctx_t ctx)
{
    if (A->length == 0)
    {
        return 0;
    }
    else if (A->bits <= FLINT_BITS)
    {
        return _nmod_mpoly_evaluate_all_ui_sp(A, vals, ctx);
    }
    else
    {
        return _nmod_mpoly_evaluate_all_ui_mp(A, vals, ctx);
    }
}

