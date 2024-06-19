/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"

ulong fmpz_mpoly_evaluate_all_nmod(
    const fmpz_mpoly_t A,
    const ulong * alphas,
    const fmpz_mpoly_ctx_t ctx,
    nmod_t fpctx)
{
    ulong eval, * t;
    TMP_INIT;

    TMP_START;

    t = TMP_ARRAY_ALLOC(A->length, ulong);
    _fmpz_vec_get_nmod_vec(t, A->coeffs, A->length, fpctx);
    eval = _nmod_mpoly_eval_all_ui(t, A->exps, A->length, A->bits,
                                                    alphas, ctx->minfo, fpctx);
    TMP_END;

    return eval;
}
