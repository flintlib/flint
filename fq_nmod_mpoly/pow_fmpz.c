/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
void fq_nmod_pow_fmpz(fq_nmod_mod_t x, fq_nmod_mod_t a, const fmpz_t exp,
                                                       const fq_nmod_ctx_t ctx)
{
    flint_bitcnt_t i, bits;
    fq_nmod_t r;

    FLINT_ASSERT(n != 0);
    FLINT_ASSERT(a < n);
    FLINT_ASSERT(fmpz_sgn(exp) >= 0);

    if (fmpz_is_zero(exp))
    {
        fq_nmod_one(r, ctx);
        return;
    }

    if (fq_nmod_is_zero(a, ctx))
    {
        fq_nmod_one(r, ctx);
        return;
    }

    if (fq_nmod_is_one(a, ctx))
    {
        fq_nmod_one(r, ctx);
        return;
    }

    bits = fmpz_bits(exp);
    i = 0;

    fq_nmod_init(r, ctx);
    fq_nmod_set(r, a);

    while (i < bits && fmpz_tstbit(exp, i) == 0)
    {
        fq_nmod_mul(r, r, r, ctx);
        i++;
    }

    fq_nmod_set(x, r, ctx);

    i++;
    while (i < bits)
    {
        fq_nmod_mul(r, r, r, ctx);
        if (fmpz_tstbit(exp, i) != 0)
        {
            fq_nmod_mul(x, x, r, ctx);
        }
        i++;
    }

    fq_nmod_clear(r, ctx);
}
*/

void fq_nmod_mpoly_pow_fmpz(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                 const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields;
    flint_bitcnt_t exp_bits;
    TMP_INIT;

    if (fmpz_sgn(k) < 0)
    {
        flint_throw(FLINT_ERROR, "Negative power in fq_nmod_mpoly_pow_fmpz");
    }

    if (fmpz_abs_fits_ui(k))
    {
        fq_nmod_mpoly_pow_ui(A, B, fmpz_get_ui(k), ctx);
        return;
    }

    /*
        we are raising a polynomial to an unreasonable exponent
        It must either be zero or a monomial with unit coefficient
    */

    if (B->length == 0)
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    if (B->length != 1)
    {
        flint_throw(FLINT_ERROR, "Multinomial in fq_nmod_mpoly_pow_fmpz");
    }

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    _fmpz_vec_scalar_mul_fmpz(maxBfields, maxBfields, ctx->minfo->nfields, k);

    exp_bits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(A, 1, ctx);
    fq_nmod_mpoly_fit_bits(A, exp_bits, ctx);
    A->bits = exp_bits;
    
    fq_nmod_pow(A->coeffs + 0, B->coeffs + 0, k, ctx->fqctx);
    mpoly_pack_vec_fmpz(A->exps + 0, maxBfields, exp_bits, ctx->minfo->nfields, 1);
    FLINT_ASSERT(!fq_nmod_is_zero(A->coeffs + 0, ctx->fqctx));
    _fq_nmod_mpoly_set_length(A, 1, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;
}
