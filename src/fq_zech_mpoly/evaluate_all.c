/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void _fq_zech_mpoly_eval_all_fq_zech(
    fq_zech_t eval,
    const fq_zech_struct * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    fq_zech_struct * const * alphas,
    const mpoly_ctx_t mctx,
    const fq_zech_ctx_t fqctx)
{
    slong i, j;
    slong nvars = mctx->nvars;
    ulong mask = (Abits <= FLINT_BITS) ? (-UWORD(1)) >> (FLINT_BITS - Abits) : 0;
    slong N = mpoly_words_per_exp(Abits, mctx);
    ulong varexp_sp;
    fmpz_t varexp_mp;
    slong * offsets, * shifts;
    fq_zech_t t, p;
    TMP_INIT;

    TMP_START;

    fmpz_init(varexp_mp);
    fq_zech_init(t, fqctx);
    fq_zech_init(p, fqctx);

    offsets = (slong *) TMP_ALLOC(2*nvars*sizeof(slong));
    shifts = offsets + nvars;
    for (j = 0; j < nvars; j++)
    {
        if (Abits <= FLINT_BITS)
            mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, Abits, mctx);
        else
            offsets[j] = mpoly_gen_offset_mp(j, Abits, mctx);
    }

    fq_zech_zero(eval, fqctx);

    for (i = 0; i < Alen; i++)
    {
        fq_zech_set(t, Acoeffs + i, fqctx);
        if (Abits <= FLINT_BITS)
        {
            for (j = 0; j < nvars; j++)
            {
                varexp_sp = ((Aexps + N*i)[offsets[j]]>>shifts[j])&mask;
                fq_zech_pow_ui(p, alphas[j], varexp_sp, fqctx);
                fq_zech_mul(t, t, p, fqctx);
            }
        }
        else
        {
            for (j = 0; j < nvars; j++)
            {
                fmpz_set_ui_array(varexp_mp, Aexps + N*i + offsets[j], Abits/FLINT_BITS);
                fq_zech_pow(p, alphas[j], varexp_mp, fqctx);
                fq_zech_mul(t, t, p, fqctx);
            }
        }

        fq_zech_add(eval, eval, t, fqctx);
    }

    fmpz_clear(varexp_mp);
    fq_zech_clear(t, fqctx);
    fq_zech_clear(p, fqctx);

    TMP_END;
}


void fq_zech_mpoly_evaluate_all_ui(fq_zech_t eval, const fq_zech_mpoly_t A,
                  fq_zech_struct * const * vals, const fq_zech_mpoly_ctx_t ctx)
{
    if (fq_zech_mpoly_is_zero(A, ctx))
    {
        fq_zech_zero(eval, ctx->fqctx);
        return;
    }

    _fq_zech_mpoly_eval_all_fq_zech(eval, A->coeffs, A->exps, A->length,
                                        A->bits, vals, ctx->minfo, ctx->fqctx);
}

