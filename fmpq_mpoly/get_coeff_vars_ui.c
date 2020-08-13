/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_get_coeff_vars_ui(fmpq_mpoly_t C, const fmpq_mpoly_t A,
                         const slong * vars, const ulong * exps, slong length,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_get_coeff_vars_ui(C->zpoly, A->zpoly, vars, exps, length, ctx->zctx);
    fmpq_set(C->content, A->content);
    fmpq_mpoly_reduce(C, ctx);
}
