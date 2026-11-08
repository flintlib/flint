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
#include "mp_real.h"
#include "mp_real/impl.h"
#include "fmpq.h"

/* the limb-radix kernel must give the exact rational
   T / (Q B^QE) = sum_{k=1}^N (x B^-D)^k / k!, here by Horner's rule
   in fmpq arithmetic */
static void
_ref_sum(fmpq_t s, const fmpz_t x, slong D, slong N)
{
    fmpq_t t, xr;
    slong k;

    fmpq_init(t);
    fmpq_init(xr);
    fmpz_set(fmpq_numref(xr), x);
    fmpz_one(fmpq_denref(xr));
    fmpz_mul_2exp(fmpq_denref(xr), fmpq_denref(xr), (ulong) (FLINT_BITS * D));
    fmpq_canonicalise(xr);

    /* s = x/1 (1 + x/2 (1 + x/3 (... (1 + x/N)))) */
    fmpq_zero(s);
    for (k = N; k >= 1; k--)
    {
        fmpq_add_ui(s, s, 1);
        fmpq_mul(s, s, xr);
        fmpz_mul_ui(fmpq_denref(s), fmpq_denref(s), (ulong) k);
        fmpq_canonicalise(s);
    }

    fmpq_clear(t);
    fmpq_clear(xr);
}

TEST_FUNCTION_START(mp_real_exp_sum_bs, state)
{
    slong iter;

    for (iter = 0; iter < 200; iter++)
    {
        slong D = 1 + n_randint(state, 8);
        slong xn = 1 + n_randint(state, D);
        slong N = 1 + n_randint(state, 60);
        slong tn, qn, QE, i;
        nn_ptr xp, T, Q;
        fmpz_t xf, num, den;
        fmpq_t qa, qb;

        /* large trees reach the threaded subtrees and merges */
        if (iter % 20 == 0)
        {
            D = 1 + n_randint(state, 40);
            xn = 1 + n_randint(state, D);
            N = 64 + n_randint(state, 400);
        }
        flint_set_num_threads(1 + n_randint(state, 4));

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

        _mp_real_exp_sum_bs_powtab(T, &tn, Q, &qn, &QE, xp, xn, D, N);

        fmpz_init(xf);
        fmpz_init(num); fmpz_init(den);
        fmpq_init(qa); fmpq_init(qb);

        fmpz_set_ui_array(xf, xp, xn);
        _ref_sum(qb, xf, D, N);

        fmpz_set_ui_array(num, T, FLINT_MAX(tn, 1));
        fmpz_set_ui_array(den, Q, FLINT_MAX(qn, 1));
        fmpz_mul_2exp(den, den, (ulong) (FLINT_BITS * QE));
        fmpq_set_fmpz_frac(qa, num, den);

        if (!fmpq_equal(qa, qb))
        {
            flint_printf("FAIL: value mismatch D=%wd xn=%wd N=%wd\n",
                D, xn, N);
            flint_abort();
        }

        fmpz_clear(xf);
        fmpz_clear(num); fmpz_clear(den);
        fmpq_clear(qa); fmpq_clear(qb);
        flint_free(xp); flint_free(T); flint_free(Q);
    }

    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
