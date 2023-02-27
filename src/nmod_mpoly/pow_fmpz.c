/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int nmod_mpoly_pow_fmpz(nmod_mpoly_t A, const nmod_mpoly_t B,
                                  const fmpz_t k, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields;
    flint_bitcnt_t exp_bits;
    TMP_INIT;

    if (fmpz_sgn(k) < 0)
        flint_throw(FLINT_ERROR, "nmod_mpoly_pow_fmpz: power is negative");

    if (fmpz_fits_si(k))
        return nmod_mpoly_pow_ui(A, B, fmpz_get_ui(k), ctx);

    /*
        we are raising a polynomial to an unreasonable exponent
        It must either be zero or a monomial with unit coefficient
    */

    if (B->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return 1;
    }

    if (B->length != 1)
        return 0;

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    _fmpz_vec_scalar_mul_fmpz(maxBfields, maxBfields, ctx->minfo->nfields, k);

    exp_bits = 1 + _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    nmod_mpoly_fit_length_reset_bits(A, 1, exp_bits, ctx);
    
    A->coeffs[0] = nmod_pow_fmpz(B->coeffs[0], k, ctx->mod);
    mpoly_pack_vec_fmpz(A->exps + 0, maxBfields, exp_bits, ctx->minfo->nfields, 1);
    _nmod_mpoly_set_length(A, A->coeffs[0] != 0, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;

    return 1;
}
