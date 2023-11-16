/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"
#include "fexpr.h"

void
fexpr_set_fmpz_mpoly_q(fexpr_t res, const fmpz_mpoly_q_t frac, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(frac), ctx))
    {
        fexpr_set_fmpz_mpoly(res, fmpz_mpoly_q_numref(frac), vars, ctx);
    }
    else
    {
        fexpr_t p, q;

        fexpr_init(p);
        fexpr_init(q);

        fexpr_set_fmpz_mpoly(p, fmpz_mpoly_q_numref(frac), vars, ctx);
        fexpr_set_fmpz_mpoly(q, fmpz_mpoly_q_denref(frac), vars, ctx);
        fexpr_div(res, p, q);

        fexpr_clear(p);
        fexpr_clear(q);

    }
}
