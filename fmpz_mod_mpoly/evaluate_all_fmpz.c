/*
    Copyright (C) 2020, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void _fmpz_mod_mpoly_eval_all_fmpz_mod(
    fmpz_t eval,
    const fmpz * Acoeffs,   /* need not be reduced mod fctx */
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const fmpz * alphas,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fctx)
{
    slong i, j;
    slong nvars = mctx->nvars;
    ulong mask = (Abits <= FLINT_BITS) ? (-UWORD(1)) >> (FLINT_BITS - Abits) : 0;
    slong N = mpoly_words_per_exp(Abits, mctx);
    ulong varexp_sp;
    fmpz_t varexp_mp, m, p;
    slong * offsets, * shifts;
    TMP_INIT;

    TMP_START;

    fmpz_init(varexp_mp);
    fmpz_init(m);
    fmpz_init(p);

    offsets = (slong *) TMP_ALLOC(2*nvars*sizeof(slong));
    shifts = offsets + nvars;
    for (j = 0; j < nvars; j++)
    {
        if (Abits <= FLINT_BITS)
            mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, Abits, mctx);
        else
            offsets[j] = mpoly_gen_offset_mp(j, Abits, mctx);
    }

    fmpz_zero(eval);
    for (i = 0; i < Alen; i++)
    {
        fmpz_one(m);
        if (Abits <= FLINT_BITS)
        {
            for (j = 0; j < nvars; j++)
            {
                varexp_sp = ((Aexps + N*i)[offsets[j]]>>shifts[j])&mask;
                fmpz_mod_pow_ui(p, alphas + j, varexp_sp, fctx);
                fmpz_mod_mul(m, m, p, fctx);
            }
        }
        else
        {
            for (j = 0; j < nvars; j++)
            {
                fmpz_set_ui_array(varexp_mp, Aexps + N*i + offsets[j], Abits/FLINT_BITS);
                fmpz_mod_pow_fmpz(p, alphas + j, varexp_mp, fctx);
                fmpz_mod_mul(m, m, p, fctx);
            }
        }

        fmpz_addmul(eval, Acoeffs + i, m);
    }

    fmpz_clear(varexp_mp);
    fmpz_clear(m);
    fmpz_clear(p);

    TMP_END;

    fmpz_mod_set_fmpz(eval, eval, fctx);
}


void fmpz_mod_mpoly_evaluate_all_fmpz(
    fmpz_t eval,
    const fmpz_mod_mpoly_t A,
    fmpz * const * vals,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * t;
    TMP_INIT;

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_zero(eval);
        return;
    }

    TMP_START;
    t = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_init(t + i);
        fmpz_mod_set_fmpz(t + i, vals[i], ctx->ffinfo);
    }

    _fmpz_mod_mpoly_eval_all_fmpz_mod(eval, A->coeffs, A->exps, A->length,
                                          A->bits, t, ctx->minfo, ctx->ffinfo);

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(t + i);

    TMP_END;
}
