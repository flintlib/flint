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

/* the joint truncated splitting kernel: for random D, xn <= D, N
   and lmax, the returned integers satisfy

       A B^ae / (Q B^QE)                     ~ 1 - cos(x B^-D),
       (x B^-D) (Q B^QE - B B^be) / (Q B^QE) ~ sin(x B^-D)

   within the truncated tree's working accuracy (relative
   2^(-64 (lmax - 1) + 40) covers the one-sided drop class with
   margin) plus the series tail beyond N terms */

TEST_FUNCTION_START(fixed_sin_cos_sum_bs, state)
{
    slong iter;

    for (iter = 0; iter < 50 + 100 * flint_test_multiplier();
        iter++)
    {
        slong D = 1 + n_randint(state, 60);
        slong xn = 1 + n_randint(state, D);
        slong xbits = 64 * (xn - 1) + 1
            + (slong) n_randint(state, 64);
        slong N = 1 + n_randint(state, 60);
        slong lmax = 4 + n_randint(state, 40);
        nn_ptr x, A, B, Q;
        slong an, bn, qn, ae, be, qb2, QE;
        fmpz_t fx, fa, fb, fq;
        arb_t va, vb, ra, rb, t, tol, tail;
        slong prec = 128 * D + 64 * lmax + 64 * N + 512;

        x = flint_malloc((xn + 1) * sizeof(ulong));
        flint_mpn_urandomb(x, state, (flint_bitcnt_t) xbits);
        while (xn > 1 && x[xn - 1] == 0)
            xn--;
        if (x[xn - 1] == 0)
            x[xn - 1] = 1;

        qb2 = (N * 2 * FLINT_BIT_COUNT(2 * (ulong) N + 1))
            / FLINT_BITS + 3;
        A = flint_malloc((lmax + qb2 + 8) * sizeof(ulong));
        B = flint_malloc((lmax + qb2 + 8) * sizeof(ulong));
        Q = flint_malloc(qb2 * sizeof(ulong));

        _fixed_sin_cos_sum_bs_powtab(A, &an, &ae, B, &bn, &be,
            Q, &qn, &QE, x, xn, D, N, lmax);

        fmpz_init(fx); fmpz_init(fa); fmpz_init(fb); fmpz_init(fq);
        arb_init(va); arb_init(vb); arb_init(ra); arb_init(rb);
        arb_init(t); arb_init(tol); arb_init(tail);

        fmpz_set_ui_array(fx, x, xn);
        fmpz_set_ui_array(fa, A, FLINT_MAX(an, 1));
        fmpz_set_ui_array(fb, B, FLINT_MAX(bn, 1));
        fmpz_set_ui_array(fq, Q, qn);

        /* va = A B^ae / (Q B^QE) ~ 1 - cos(x B^-D) */
        arb_set_fmpz(va, fa);
        arb_mul_2exp_si(va, va, 64 * ae);
        arb_set_fmpz(t, fq);
        arb_mul_2exp_si(t, t, 64 * QE);
        arb_div(va, va, t, prec);
        /* vb = (x B^-D)(Q B^QE - B B^be)/(Q B^QE) ~ sin */
        arb_set_fmpz(vb, fb);
        arb_mul_2exp_si(vb, vb, 64 * be);
        arb_sub(vb, t, vb, prec);
        arb_div(vb, vb, t, prec);
        arb_set_fmpz(t, fx);
        arb_mul_2exp_si(t, t, -64 * D);
        arb_mul(vb, vb, t, prec);

        arb_sin_cos(rb, ra, t, prec);
        arb_sub_ui(ra, ra, 1, prec);
        arb_neg(ra, ra);

        arb_sub(ra, ra, va, prec);
        arb_abs(ra, ra);
        arb_sub(rb, rb, vb, prec);
        arb_abs(rb, rb);

        arb_abs(tol, va);
        arb_mul_2exp_si(tol, tol, -64 * (lmax - 1) + 40);
        arb_set_fmpz(tail, fx);
        arb_mul_2exp_si(tail, tail, -64 * D);
        arb_pow_ui(tail, tail, 2 * (ulong) N + 2, 128);
        arb_add(tol, tol, tail, 128);
        if (!arb_le(ra, tol))
            TEST_FUNCTION_FAIL("A track: D = %wd, xbits = %wd, "
                "N = %wd, lmax = %wd\n", D, xbits, N, lmax);

        arb_abs(tol, vb);
        arb_mul_2exp_si(tol, tol, -64 * (lmax - 1) + 40);
        arb_set_fmpz(tail, fx);
        arb_mul_2exp_si(tail, tail, -64 * D);
        arb_pow_ui(tail, tail, 2 * (ulong) N + 1, 128);
        arb_add(tol, tol, tail, 128);
        if (!arb_le(rb, tol))
            TEST_FUNCTION_FAIL("B track: D = %wd, xbits = %wd, "
                "N = %wd, lmax = %wd\n", D, xbits, N, lmax);

        flint_free(x); flint_free(A); flint_free(B); flint_free(Q);
        fmpz_clear(fx); fmpz_clear(fa); fmpz_clear(fb);
        fmpz_clear(fq);
        arb_clear(va); arb_clear(vb); arb_clear(ra); arb_clear(rb);
        arb_clear(t); arb_clear(tol); arb_clear(tail);
    }

    TEST_FUNCTION_END(state);
}
