/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"

void
fmpz_mpoly_q_factored_randtest(fmpz_mpoly_q_factored_t res, flint_rand_t state,
    slong length, ulong coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_struct * num = fmpz_mpoly_q_factored_numref(res);
    fmpz_mpoly_factor_struct * den = fmpz_mpoly_q_factored_denref(res);
    fmpz_mpoly_t t;
    slong i, nfac;

    length = FLINT_MAX(length, 1);
    exp_bound = FLINT_MAX(exp_bound, 1);
    coeff_bits = FLINT_MAX(coeff_bits, 1);

    fmpz_mpoly_randtest_bound(num, state, length, coeff_bits, exp_bound, ctx);

    fmpz_mpoly_factor_one(den, ctx);
    fmpz_one(den->constant_den);

    if (n_randint(state, 2))
    {
        fmpz_randtest_not_zero(den->constant, state, coeff_bits);
        fmpz_abs(den->constant, den->constant);
    }

    fmpz_mpoly_init(t, ctx);

    nfac = n_randint(state, 5);
    for (i = 0; i < nfac; i++)
    {
        do {
            fmpz_mpoly_randtest_bound(t, state, 1 + n_randint(state, length),
                1 + n_randint(state, coeff_bits), exp_bound, ctx);
        } while (fmpz_mpoly_is_zero(t, ctx));

        fmpz_mpoly_factor_append_ui(den, t, 1 + n_randint(state, 3), ctx);

        /* sometimes introduce common factors to test cancellation */
        if (n_randint(state, 4) == 0)
            fmpz_mpoly_mul(num, num, t, ctx);
    }

    fmpz_mpoly_clear(t, ctx);

    fmpz_mpoly_q_factored_canonicalise(res, ctx);
}
