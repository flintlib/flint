/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "fmpz_mpoly.h"

char * fmpz_mod_mpoly_get_str_pretty(const fmpz_mod_mpoly_t A,
                               const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
{
   return _fmpz_mpoly_get_str_pretty(A->coeffs, A->exps,
                                            A->length, x, A->bits, ctx->minfo);
}
