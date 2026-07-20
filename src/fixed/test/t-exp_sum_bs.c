/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fixed.h"
#include "fmpq.h"
#include "arb.h"

/* the limb-radix kernel must agree with _arb_exp_sum_bs_powtab as
   an exact rational: T / (Q B^QE) == Ta / (Qa 2^Qexpa) */
void _arb_exp_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong N);

TEST_FUNCTION_START(fixed_exp_sum_bs, state)
{
    slong iter;

    for (iter = 0; iter < 200; iter++)
    {
        slong D = 1 + n_randint(state, 8);
        slong xn = 1 + n_randint(state, D);
        slong N = 1 + n_randint(state, 60);
        slong tn, qn, QE, i;
        nn_ptr xp, T, Q;
        fmpz_t Ta, Qa, xf, num, den;
        fmpq_t qa, qb;
        flint_bitcnt_t Qexpa;

        xp = flint_malloc(xn * sizeof(ulong));
        flint_mpn_urandomb(xp, state, FLINT_BITS * xn);
        xp[xn - 1] |= (UWORD(1) << n_randint(state, FLINT_BITS));
        while (xn > 1 && xp[0] == 0)    /* the driver strips */
        {
            for (i = 0; i < xn - 1; i++)
                xp[i] = xp[i + 1];
            xn--;
        }

        T = flint_malloc((N * (D + 2) + 4) * sizeof(ulong));
        Q = flint_malloc(((N * FLINT_BIT_COUNT((ulong) N + 1))
            / FLINT_BITS + 3) * sizeof(ulong));

        _fixed_exp_sum_bs_powtab(T, &tn, Q, &qn, &QE, xp, xn, D, N);

        fmpz_init(Ta); fmpz_init(Qa); fmpz_init(xf);
        fmpz_init(num); fmpz_init(den);
        fmpq_init(qa); fmpq_init(qb);

        fmpz_set_ui_array(xf, xp, xn);
        _arb_exp_sum_bs_powtab(Ta, Qa, &Qexpa, xf,
            (flint_bitcnt_t) (FLINT_BITS * D), N);

        fmpz_set_ui_array(num, T, FLINT_MAX(tn, 1));
        fmpz_set_ui_array(den, Q, FLINT_MAX(qn, 1));
        fmpz_mul_2exp(den, den, (ulong) (FLINT_BITS * QE));
        fmpq_set_fmpz_frac(qa, num, den);

        fmpz_mul_2exp(Qa, Qa, Qexpa);
        fmpq_set_fmpz_frac(qb, Ta, Qa);

        if (!fmpq_equal(qa, qb))
        {
            flint_printf("FAIL: value mismatch D=%wd xn=%wd N=%wd\n",
                D, xn, N);
            flint_abort();
        }

        fmpz_clear(Ta); fmpz_clear(Qa); fmpz_clear(xf);
        fmpz_clear(num); fmpz_clear(den);
        fmpq_clear(qa); fmpq_clear(qb);
        flint_free(xp); flint_free(T); flint_free(Q);
    }

    TEST_FUNCTION_END(state);
}
