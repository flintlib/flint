/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_q.h"

void
fmpz_mod_mpoly_q_init(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_init(fmpz_mod_mpoly_q_numref(res), ctx);
    fmpz_mod_mpoly_init(fmpz_mod_mpoly_q_denref(res), ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}
