/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "fmpz_mod_mpoly_q.h"

void
fmpz_mod_mpoly_q_set(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (res != x)
    {
        fmpz_mod_mpoly_set(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(x), ctx);
        fmpz_mod_mpoly_set(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(x), ctx);
        fmpz_mod_mpoly_q_canonicalise(res,ctx);
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

int
fmpz_mod_mpoly_q_set_fmpq(fmpz_mod_mpoly_q_t res, const fmpq_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_is_one(fmpq_denref(x)))
    {
        fmpz_mod_mpoly_q_set_fmpz(res, fmpq_numref(x), ctx);
        return 1;
    }
    else
    {
        fmpz_t t;
        int invertible;

        fmpz_init(t);

        fmpz_mod_set_fmpz(t, fmpq_denref(x), ctx->ffinfo);
        invertible = !fmpz_is_zero(t);

        if (invertible)
        {
            fmpz_mod_inv(t, t, ctx->ffinfo);
            fmpz_mod_mul_fmpz(t, t, fmpq_numref(x), ctx->ffinfo);
            fmpz_mod_mpoly_set_fmpz_mod(fmpz_mod_mpoly_q_numref(res), t, ctx);
            fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
        }

        fmpz_clear(t);

        return invertible;
    }
}

