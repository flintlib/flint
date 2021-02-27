/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


int fmpz_mod_mpoly_pow_ui(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    ulong k,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, exp_bits;
    fmpz * maxBfields;
    TMP_INIT;

    if (k == 0)
    {
        fmpz_mod_mpoly_one(A, ctx);
        return 1;
    }

    if (B->length == 0)
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return 1;
    }

    if (k == 1)
    {
        fmpz_mod_mpoly_set(A, B, ctx);
        return 1;
    }

    if (k == 2)
    {
        fmpz_mod_mpoly_mul(A, B, B, ctx);
        return 1;
    }

    if (B->length != 1)
    {
        fmpz_mod_mpoly_pow_rmul(A, B, k, ctx);
        return 1;
    }

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    _fmpz_vec_scalar_mul_ui(maxBfields, maxBfields, ctx->minfo->nfields, k);

    exp_bits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    exp_bits = mpoly_fix_bits(exp_bits + 1, ctx->minfo);
    fmpz_mod_mpoly_fit_length_reset_bits(A, 1, exp_bits, ctx);

    fmpz_mod_pow_ui(A->coeffs + 0, B->coeffs + 0, k, ctx->ffinfo);
    mpoly_pack_vec_fmpz(A->exps + 0, maxBfields, exp_bits, ctx->minfo->nfields, 1);
    _fmpz_mod_mpoly_set_length(A, !fmpz_is_zero(A->coeffs + 0), ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;

    return 1;
}
