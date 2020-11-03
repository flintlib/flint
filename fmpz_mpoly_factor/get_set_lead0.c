/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void _fmpz_mpoly_get_lead0(
    fmpz_mpoly_t c,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpolyl_lead_coeff(c, A, 1, ctx);
}

void _fmpz_mpoly_set_lead0(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong deg;
    fmpz_mpoly_t t, g;

    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(g, ctx);

    deg = fmpz_mpoly_degree_si(B, 0, ctx);
    FLINT_ASSERT(deg >= 0);
    fmpz_mpoly_gen(g, 0, ctx);
    fmpz_mpoly_pow_ui(g, g, deg, ctx);
    _fmpz_mpoly_get_lead0(t, B, ctx);
    fmpz_mpoly_sub(t, c, t, ctx);
    fmpz_mpoly_mul(t, t, g, ctx);
    fmpz_mpoly_add(A, B, t, ctx);

    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(g, ctx);
}
