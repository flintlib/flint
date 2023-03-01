/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

void
fmpz_mpoly_q_init(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_init(fmpz_mpoly_q_numref(res), ctx);
    fmpz_mpoly_init(fmpz_mpoly_q_denref(res), ctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), ctx);
}

