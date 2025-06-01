/*
    Copyright (C) 2020 Fredrik Johansson

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

void
fmpz_mod_mpoly_q_set_fmpq(fmpz_mod_mpoly_q_t res, const fmpq_t x, const fmpz_mod_mpoly_ctx_t ctx)
{   

    fmpz_t n, d;

    fmpz_init(n);
    fmpz_init(d);

    fmpz_mod_set_fmpz(n, fmpq_numref(x), ctx->ffinfo);
    fmpz_mod_set_fmpz(d, fmpq_denref(x), ctx->ffinfo);

    if (fmpz_is_zero(d))
    {
        fmpz_zero(d);
    }
    else
    {
        fmpz_mod_inv(d, d, ctx->ffinfo);
    }       

    fmpz_mod_mul_fmpz(n, n, d, ctx->ffinfo);
    fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_q_numref(res), n, ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);

    fmpz_clear(n);
    fmpz_clear(d);

}
