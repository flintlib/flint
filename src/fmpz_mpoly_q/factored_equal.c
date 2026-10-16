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

int
fmpz_mpoly_q_factored_equal(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)
{
    const fmpz_mpoly_factor_struct * xd = fmpz_mpoly_q_factored_denref(x);
    const fmpz_mpoly_factor_struct * yd = fmpz_mpoly_q_factored_denref(y);
    slong i;

    if (x == y)
        return 1;

    if (xd->num != yd->num || !fmpz_equal(xd->constant, yd->constant))
        return 0;

    for (i = 0; i < xd->num; i++)
        if (!fmpz_equal(xd->exp + i, yd->exp + i))
            return 0;

    if (!fmpz_mpoly_equal(fmpz_mpoly_q_factored_numref(x), fmpz_mpoly_q_factored_numref(y), ctx))
        return 0;

    for (i = 0; i < xd->num; i++)
        if (!fmpz_mpoly_equal(xd->poly + i, yd->poly + i, ctx))
            return 0;

    return 1;
}
