/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_q.h"

void
fmpz_mod_mpoly_q_set(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (res != x)
    {
        fmpz_mod_mpoly_set(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(x), ctx);
        fmpz_mod_mpoly_set(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(x), ctx);
    }
}

void
fmpz_mod_mpoly_q_set_si(fmpz_mod_mpoly_q_t res, slong x, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_set_si(fmpz_mod_mpoly_q_numref(res), x, ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}

void
fmpz_mod_mpoly_q_set_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_q_numref(res), x, ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}

void
fmpz_mod_mpoly_q_set_fmpq(fmpz_mod_mpoly_q_t res, const fmpq_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_q_numref(res), fmpq_numref(x), ctx);
    fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_q_denref(res), fmpq_denref(x), ctx);
}
