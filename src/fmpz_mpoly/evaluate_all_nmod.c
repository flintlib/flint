/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/


#include "fmpz_mpoly.h"

mp_limb_t fmpz_mpoly_evaluate_all_nmod(
    const fmpz_mpoly_t A,
    const mp_limb_t * alphas,
    const fmpz_mpoly_ctx_t ctx,
    nmod_t fpctx)
{
    mp_limb_t eval, * t;
    TMP_INIT;

    TMP_START;

    t = TMP_ARRAY_ALLOC(A->length, mp_limb_t);
    _fmpz_vec_get_nmod_vec(t, A->coeffs, A->length, fpctx);
    eval = _nmod_mpoly_eval_all_ui(t, A->exps, A->length, A->bits,
                                                    alphas, ctx->minfo, fpctx);
    TMP_END;

    return eval;
}

