/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"

void
fmpz_mpoly_q_factored_init(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_init(fmpz_mpoly_q_factored_numref(res), ctx);
    fmpz_mpoly_factor_init(fmpz_mpoly_q_factored_denref(res), ctx);
}

void
fmpz_mpoly_q_factored_clear(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_clear(fmpz_mpoly_q_factored_numref(res), ctx);
    fmpz_clear(res->den.constant_den);
    fmpz_mpoly_factor_clear(fmpz_mpoly_q_factored_denref(res), ctx);
}

void
fmpz_mpoly_q_factored_set(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
{
    if (res != x)
    {
        fmpz_mpoly_set(fmpz_mpoly_q_factored_numref(res), fmpz_mpoly_q_factored_numref(x), ctx);
        fmpz_mpoly_factor_set(fmpz_mpoly_q_factored_denref(res), fmpz_mpoly_q_factored_denref(x), ctx);
    }
}

void
fmpz_mpoly_q_factored_zero(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_zero(fmpz_mpoly_q_factored_numref(res), ctx);
    fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
}

void
fmpz_mpoly_q_factored_one(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_one(fmpz_mpoly_q_factored_numref(res), ctx);
    fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
}

void
fmpz_mpoly_q_factored_gen(fmpz_mpoly_q_factored_t res, slong i, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_gen(fmpz_mpoly_q_factored_numref(res), i, ctx);
    fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
}

void
fmpz_mpoly_q_factored_set_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set(fmpz_mpoly_q_factored_numref(res), x, ctx);
    fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
}

void
fmpz_mpoly_q_factored_set_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_factored_numref(res), x, ctx);
    fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
}

void
fmpz_mpoly_q_factored_set_si(fmpz_mpoly_q_factored_t res, slong x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_si(fmpz_mpoly_q_factored_numref(res), x, ctx);
    fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
}

void
fmpz_mpoly_q_factored_set_fmpq(fmpz_mpoly_q_factored_t res, const fmpq_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_factored_numref(res), fmpq_numref(x), ctx);
    fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
    if (!fmpz_is_zero(fmpq_numref(x)))
        fmpz_set(res->den.constant, fmpq_denref(x));
}

/* Converting from fmpz_mpoly_q_t requires factoring the denominator.
   Since the input is canonical, no cancellation is necessary. */
void
fmpz_mpoly_q_factored_set_fmpz_mpoly_q(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_factor_t f;
    fmpz_t s;

    fmpz_mpoly_factor_init(f, ctx);
    fmpz_init(s);

    _fmpz_mpoly_q_factored_factor_num(s, f, fmpz_mpoly_q_denref(x), ctx);

    if (fmpz_sgn(s) < 0)
    {
        fmpz_neg(s, s);
        fmpz_mpoly_neg(fmpz_mpoly_q_factored_numref(res), fmpz_mpoly_q_numref(x), ctx);
    }
    else
    {
        fmpz_mpoly_set(fmpz_mpoly_q_factored_numref(res), fmpz_mpoly_q_numref(x), ctx);
    }

    fmpz_swap(f->constant, s);
    fmpz_mpoly_factor_swap(fmpz_mpoly_q_factored_denref(res), f, ctx);

    fmpz_mpoly_factor_clear(f, ctx);
    fmpz_clear(s);
}

/* Since the input is canonical, the expanded denominator is automatically
   coprime to the numerator and has positive leading coefficient. */
void
fmpz_mpoly_q_factored_get_fmpz_mpoly_q(fmpz_mpoly_q_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_factored_numref(x), ctx);
    if (!fmpz_mpoly_factor_expand(fmpz_mpoly_q_denref(res), fmpz_mpoly_q_factored_denref(x), ctx))
        flint_throw(FLINT_ERROR, "fmpz_mpoly_q_factored_get_fmpz_mpoly_q: expansion failed\n");
}
