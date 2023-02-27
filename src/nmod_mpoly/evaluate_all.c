/*
    Copyright (C) 2018,2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

mp_limb_t _nmod_mpoly_eval_all_ui(
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const mp_limb_t * alphas,
    const mpoly_ctx_t mctx,
    nmod_t mod)
{
    slong i, j;
    slong nvars = mctx->nvars;
    ulong mask = (Abits <= FLINT_BITS) ? (-UWORD(1)) >> (FLINT_BITS - Abits) : 0;
    slong N = mpoly_words_per_exp(Abits, mctx);
    ulong varexp_sp;
    fmpz_t varexp_mp;
    slong * offsets, * shifts;
    n_poly_struct * caches;
    mp_limb_t eval, t;
    TMP_INIT;

    TMP_START;

    fmpz_init(varexp_mp);

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

        t = alphas[j];
        if (t >= mod.n)
            NMOD_RED(t, t, mod);

        nmod_pow_cache_start(t, caches + 3*j + 0, caches + 3*j + 1,
                                                          caches + 3*j + 2);
    }

    eval = 0;
    for (i = 0; i < Alen; i++)
    {
        t = Acoeffs[i];
        if (Abits <= FLINT_BITS)
        {
            for (j = 0; j < nvars; j++)
            {
                varexp_sp = ((Aexps + N*i)[offsets[j]]>>shifts[j])&mask;
                t = nmod_pow_cache_mulpow_ui(t, varexp_sp, caches + 3*j + 0,
                                      caches + 3*j + 1, caches + 3*j + 2, mod);
            }
        }
        else
        {
            for (j = 0; j < nvars; j++)
            {
                fmpz_set_ui_array(varexp_mp, Aexps + N*i + offsets[j], Abits/FLINT_BITS);
                t = nmod_pow_cache_mulpow_fmpz(t, varexp_mp, caches + 3*j + 0,
                                      caches + 3*j + 1, caches + 3*j + 2, mod);
            }
        }

        eval = nmod_add(eval, t, mod);
    }

    fmpz_clear(varexp_mp);

    for (j = 0; j < 3*nvars; j++)
        n_poly_clear(caches + j);

    TMP_END;

    return eval;
}


ulong nmod_mpoly_evaluate_all_ui(const nmod_mpoly_t A,
                                const ulong * vals, const nmod_mpoly_ctx_t ctx)
{
    if (nmod_mpoly_is_zero(A, ctx))
        return 0;

    return _nmod_mpoly_eval_all_ui(A->coeffs, A->exps, A->length, A->bits,
                                                   vals, ctx->minfo, ctx->mod);
}

