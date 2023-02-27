/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_ctx_init(fmpz_mod_mpoly_ctx_t ctx, slong nvars,
                                    const ordering_t ord, const fmpz_t modulus)
{
    mpoly_ctx_init(ctx->minfo, nvars, ord);
    fmpz_mod_ctx_init(ctx->ffinfo, modulus);
}

