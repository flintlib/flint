/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_q.h"

void
fmpz_mod_mpoly_q_randtest(fmpz_mod_mpoly_q_t res, flint_rand_t state,
    slong length, slong exp_bound, const fmpz_mod_mpoly_ctx_t ctx)
{
    length = n_randint(state, length + 1);

    fmpz_mod_mpoly_randtest_bound(fmpz_mod_mpoly_q_numref(res), state, length, exp_bound, ctx);

    if (n_randint(state, 2))
    {
        fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
    }
    else
    {
        if (n_randint(state, 2))
        {
            length = 1;
            exp_bound = 1;
        }

        fmpz_mod_mpoly_randtest_bound(fmpz_mod_mpoly_q_denref(res), state, length, exp_bound, ctx);

        if (fmpz_mod_mpoly_is_zero(fmpz_mod_mpoly_q_denref(res), ctx))
            fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
    }
    
    fmpz_mod_mpoly_q_canonicalise(res, ctx);
}
