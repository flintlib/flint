/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "fixed.h"

/* the mpn port must reproduce _arb_exp_sum_bs_powtab bit for bit */

TEST_FUNCTION_START(fixed_exp_sum_bs, state)
{
    slong iter;

    for (iter = 0; iter < 40 + 40 * flint_test_multiplier(); iter++)
    {
        slong xn = 1 + n_randint(state, 6);
        flint_bitcnt_t r = FLINT_BITS * (xn + n_randint(state, 4));
        slong N = 1 + n_randint(state, 300);
        slong tn, qn, i;
        flint_bitcnt_t Qexp;
        nn_ptr xp, T, Q;
        fmpz_t x, Ta, Qa, Tm, Qm;
        flint_bitcnt_t Qexpa;

        xp = flint_malloc(xn * sizeof(ulong));
        T = flint_malloc(((N * ((slong) r + 128)) / FLINT_BITS + 4)
            * sizeof(ulong));
        Q = flint_malloc(((N * FLINT_BIT_COUNT((ulong) N + 1))
            / FLINT_BITS + 3) * sizeof(ulong));

        flint_mpn_urandomb(xp, state, FLINT_BITS * xn);
        xp[xn - 1] |= UWORD(1);   /* nonzero */
        while (xp[xn - 1] == 0)
            xn--;

        fmpz_init(x); fmpz_init(Ta); fmpz_init(Qa);
        fmpz_init(Tm); fmpz_init(Qm);
        fmpz_set_ui_array(x, xp, xn);

        _arb_exp_sum_bs_powtab(Ta, Qa, &Qexpa, x, r, N);
        _fixed_exp_sum_bs_powtab(T, &tn, Q, &qn, &Qexp, xp, xn, r, N);

        fmpz_set_ui_array(Tm, T, FLINT_MAX(tn, 1));
        fmpz_set_ui_array(Qm, Q, FLINT_MAX(qn, 1));

        if (!fmpz_equal(Tm, Ta) || !fmpz_equal(Qm, Qa)
            || Qexp != Qexpa)
            TEST_FUNCTION_FAIL("xn = %wd, r = %wd, N = %wd "
                "(Qexp %wd vs %wd)\n", xn, (slong) r, N,
                (slong) Qexp, (slong) Qexpa);

        /* lengths must be normalized */
        if ((tn > 0 && T[tn - 1] == 0) || (qn > 0 && Q[qn - 1] == 0))
            TEST_FUNCTION_FAIL("unnormalized output\n");
        (void) i;

        fmpz_clear(x); fmpz_clear(Ta); fmpz_clear(Qa);
        fmpz_clear(Tm); fmpz_clear(Qm);
        flint_free(xp); flint_free(T); flint_free(Q);
    }

    /* the terms helper agrees with the arb pair it ports, up to the
       safe floor-log continuation beyond the table */
    {
        slong r2, p2;
        for (r2 = 16; r2 <= 512; r2 *= 2)
            for (p2 = 128; p2 <= 65536; p2 *= 4)
            {
                slong Nf = _fixed_exp_bs_num_terms(
                    (flint_bitcnt_t) r2, p2);
                slong Na = _arb_exp_taylor_bound(-r2, p2) - 1;
                if (Na > 10000) while (Na % 128) Na++;
                if (Na > 1000) while (Na % 16) Na++;
                if (Na > 100) while (Na % 2) Na++;
                if (Nf < Na || Nf > Na + Na / 8 + 16)
                    TEST_FUNCTION_FAIL("terms: r = %wd, prec = %wd: "
                        "%wd vs arb %wd\n", r2, p2, Nf, Na);
            }
    }

    TEST_FUNCTION_END(state);
}
