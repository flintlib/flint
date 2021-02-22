/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "fmpz_mod_mpoly.h"

void fmpz_mpoly_evaluate_all_fmpz_mod(
    fmpz_t ev,
    const fmpz_mpoly_t A,
    const fmpz * alphas,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    _fmpz_mod_mpoly_eval_all_fmpz_mod(ev, A->coeffs, A->exps, A->length,
                 A->bits, alphas, ctx->minfo, fpctx);
}

