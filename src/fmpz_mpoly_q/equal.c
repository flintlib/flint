/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

int
fmpz_mpoly_q_equal(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_equal(fmpz_mpoly_q_numref(x), fmpz_mpoly_q_numref(y), ctx) &&
            fmpz_mpoly_equal(fmpz_mpoly_q_denref(x), fmpz_mpoly_q_denref(y), ctx);
}

