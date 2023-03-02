/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_divrem_ideal(
    fmpz_mod_mpoly_struct ** Q,
    fmpz_mod_mpoly_t R,
    const fmpz_mod_mpoly_t A,
    fmpz_mod_mpoly_struct * const * B,
    slong len,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_divrem_ideal_monagan_pearce(Q, R, A, B, len, ctx);
}
