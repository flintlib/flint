/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

void
fmpz_mpoly_q_set(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
{
    if (res != x)
    {
        fmpz_mpoly_set(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(x), ctx);
        fmpz_mpoly_set(fmpz_mpoly_q_denref(res), fmpz_mpoly_q_denref(x), ctx);
    }
}

void
fmpz_mpoly_q_set_si(fmpz_mpoly_q_t res, slong x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_si(fmpz_mpoly_q_numref(res), x, ctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), ctx);
}

void
fmpz_mpoly_q_set_fmpz(fmpz_mpoly_q_t res, const fmpz_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_numref(res), x, ctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), ctx);
}

void
fmpz_mpoly_q_set_fmpq(fmpz_mpoly_q_t res, const fmpq_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_numref(res), fmpq_numref(x), ctx);
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), fmpq_denref(x), ctx);
}
