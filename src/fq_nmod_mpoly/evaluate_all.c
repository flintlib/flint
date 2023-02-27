/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void _fq_nmod_mpoly_eval_all_fq_nmod(
    fq_nmod_t eval,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    fq_nmod_struct * const * alphas,
    const mpoly_ctx_t mctx,
    const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
    slong i, j;
    slong nvars = mctx->nvars;
    ulong mask = (Abits <= FLINT_BITS) ? (-UWORD(1)) >> (FLINT_BITS - Abits) : 0;
    slong N = mpoly_words_per_exp(Abits, mctx);
    ulong varexp_sp;
    fmpz_t varexp_mp;
    slong * offsets, * shifts;
    n_poly_struct * caches;
    mp_limb_t * t;
    TMP_INIT;

    TMP_START;

    fmpz_init(varexp_mp);

    t = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));
    caches = (n_poly_struct *) TMP_ALLOC(3*nvars*sizeof(n_poly_struct));
    offsets = (slong *) TMP_ALLOC(2*nvars*sizeof(slong));
    shifts = offsets + nvars;
    for (j = 0; j < nvars; j++)
    {
        if (Abits <= FLINT_BITS)
            mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, Abits, mctx);
        else
            offsets[j] = mpoly_gen_offset_mp(j, Abits, mctx);

        n_poly_init(caches + 3*j + 0);
        n_poly_init(caches + 3*j + 1);
        n_poly_init(caches + 3*j + 2);

        n_fq_pow_cache_start_fq_nmod(alphas[j], caches + 3*j + 0,
                                    caches + 3*j + 1, caches + 3*j + 2, fqctx);
    }

    nmod_poly_fit_length(eval, d);
    _n_fq_zero(eval->coeffs, d);
    for (i = 0; i < Alen; i++)
    {
        _n_fq_set(t, Acoeffs + d*i, d);
        if (Abits <= FLINT_BITS)
        {
            for (j = 0; j < nvars; j++)
            {
                varexp_sp = ((Aexps + N*i)[offsets[j]]>>shifts[j])&mask;
                n_fq_pow_cache_mulpow_ui(t, t, varexp_sp, caches + 3*j + 0,
                                    caches + 3*j + 1, caches + 3*j + 2, fqctx);
            }
        }
        else
        {
            for (j = 0; j < nvars; j++)
            {
                fmpz_set_ui_array(varexp_mp, Aexps + N*i + offsets[j], Abits/FLINT_BITS);
                n_fq_pow_cache_mulpow_fmpz(t, t, varexp_mp, caches + 3*j + 0,
                                    caches + 3*j + 1, caches + 3*j + 2, fqctx);
            }
        }

        _n_fq_add(eval->coeffs, eval->coeffs, t, d, fqctx->mod);
    }

    _nmod_poly_set_length(eval, d);
    _nmod_poly_normalise(eval);

    fmpz_clear(varexp_mp);

    for (j = 0; j < 3*nvars; j++)
        n_poly_clear(caches + j);

    TMP_END;
}


void fq_nmod_mpoly_evaluate_all_fq_nmod(fq_nmod_t ev, const fq_nmod_mpoly_t A,
                  fq_nmod_struct * const * vals, const fq_nmod_mpoly_ctx_t ctx)
{
    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_zero(ev, ctx->fqctx);
        return;
    }

    _fq_nmod_mpoly_eval_all_fq_nmod(ev, A->coeffs, A->exps, A->length, A->bits,
                                                 vals, ctx->minfo, ctx->fqctx);
}

