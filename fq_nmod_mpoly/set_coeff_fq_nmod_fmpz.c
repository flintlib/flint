/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void _fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(
    fq_nmod_mpoly_t A,
    const fq_nmod_t c,
    const fmpz * exp,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    flint_bitcnt_t exp_bits;
    slong i, N, index;
    ulong * cmpmask;
    ulong * packed_exp;
    int exists;
    TMP_INIT;

    TMP_START;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fq_nmod_mpoly_fit_length_fit_bits(A, A->length, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, A->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_set_monomial_ffmpz(packed_exp, exp, A->bits, ctx->minfo);
    exists = mpoly_monomial_exists(&index, A->exps,
                                  packed_exp, A->length, N, cmpmask);

    if (!exists)
    {
        if (!fq_nmod_is_zero(c, ctx->fqctx))
        {       
            /* make new term only if coeff is nonzero*/
            fq_nmod_mpoly_fit_length(A, A->length + 1, ctx);

            for (i = A->length; i >= index + 1; i--)
            {
                _n_fq_set(A->coeffs + d*i, A->coeffs + d*(i - 1), d);
                mpoly_monomial_set(A->exps + N*i, A->exps + N*(i - 1), N);
            }

            n_fq_set_fq_nmod(A->coeffs + d*index, c, ctx->fqctx);
            mpoly_monomial_set(A->exps + N*index, packed_exp, N);

        _fq_nmod_mpoly_set_length(A, A->length + 1, ctx);
        }
    }
    else if (fq_nmod_is_zero(c, ctx->fqctx)) /* zero coeff, remove term */
    {
        for (i = index; i < A->length - 1; i++)
        {
            _n_fq_set(A->coeffs + d*i, A->coeffs + d*(i + 1), d);
            mpoly_monomial_set(A->exps + N*i, A->exps + N*(i + 1), N);
        }

        _fq_nmod_mpoly_set_length(A, A->length - 1, ctx);
    }
    else /* term with that monomial exists, coeff is nonzero */
    {
        n_fq_set_fq_nmod(A->coeffs + d*index, c, ctx->fqctx);  
    }

    TMP_END; 
}


void fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(
    fq_nmod_mpoly_t A,
    const fq_nmod_t c,
    fmpz * const * exp,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
    {
        fmpz_init(newexp + i);
        fmpz_set(newexp + i, exp[i]);
    }

    _fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(A, c, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}
