/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mpoly_factor.h"

void fmpz_mpoly_factor_get_constant_fmpz(fmpz_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    fmpz_set(c, f->constant);
}

void fmpz_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    fmpz_set(fmpq_numref(c), f->constant);
    fmpz_set(fmpq_denref(c), f->constant_den);
}

slong fmpz_mpoly_factor_get_exp_si(fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_ASSERT(i < (ulong) f->num);
    return fmpz_get_si(f->exp + i);
}

void fmpz_mpoly_factor_set_fmpz(fmpz_mpoly_factor_t f, const fmpz_t a, const fmpz_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    f->num = 0;
    fmpz_set(f->constant, a);
}

void fmpz_mpoly_factor_zero(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    f->num = 0;
    fmpz_zero(f->constant);
}

void fmpz_mpoly_factor_one(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    f->num = 0;
    fmpz_one(f->constant);
}

void fmpz_mpoly_factor_append_fmpz_swap(fmpz_mpoly_factor_t f, fmpz_mpoly_t A, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mpoly_swap(f->poly + i, A, ctx);
    fmpz_set(f->exp + i, e);
    f->num = i + 1;
}

void fmpz_mpoly_factor_append_ui(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, ulong e, const fmpz_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mpoly_set(f->poly + i, A, ctx);
    fmpz_set_ui(f->exp + i, e);
    f->num = i + 1;
}

void fmpz_mpoly_unit_normalize(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    if (fmpz_sgn(A->coeffs + 0) < 0)
        fmpz_mpoly_neg(A, A, ctx);
}
