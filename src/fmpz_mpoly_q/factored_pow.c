/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"

/* (xn / xd)^e = xn^e / xd^e is canonical for e > 0 since powers of
   coprime polynomials are coprime. */
static void
_fmpz_mpoly_q_factored_pow_fmpz_pos(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_factor_struct * rd;
    slong i;

    fmpz_mpoly_q_factored_set(res, x, ctx);
    rd = fmpz_mpoly_q_factored_denref(res);

    if (!fmpz_mpoly_pow_fmpz(fmpz_mpoly_q_factored_numref(res), fmpz_mpoly_q_factored_numref(res), e, ctx) ||
        !fmpz_pow_fmpz(rd->constant, rd->constant, e))
    {
        flint_throw(FLINT_ERROR, "fmpz_mpoly_q_factored_pow_fmpz: exponent too large\n");
    }

    for (i = 0; i < rd->num; i++)
        fmpz_mul(rd->exp + i, rd->exp + i, e);
}

void
fmpz_mpoly_q_factored_pow_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_is_zero(e))
    {
        fmpz_mpoly_q_factored_one(res, ctx);
    }
    else if (fmpz_is_one(e))
    {
        fmpz_mpoly_q_factored_set(res, x, ctx);
    }
    else if (fmpz_mpoly_q_factored_is_zero(x, ctx))
    {
        if (fmpz_sgn(e) < 0)
            flint_throw(FLINT_DIVZERO, "fmpz_mpoly_q_factored_pow_fmpz: division by zero\n");

        fmpz_mpoly_q_factored_zero(res, ctx);
    }
    else if (fmpz_sgn(e) > 0)
    {
        _fmpz_mpoly_q_factored_pow_fmpz_pos(res, x, e, ctx);
    }
    else
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        fmpz_mpoly_q_factored_inv(res, x, ctx);
        if (!fmpz_is_one(f))
            _fmpz_mpoly_q_factored_pow_fmpz_pos(res, res, f, ctx);
        fmpz_clear(f);
    }
}

void
fmpz_mpoly_q_factored_pow_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong e, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t f;
    fmpz_init_set_si(f, e);
    fmpz_mpoly_q_factored_pow_fmpz(res, x, f, ctx);
    fmpz_clear(f);
}
