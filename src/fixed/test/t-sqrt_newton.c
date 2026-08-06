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
#include "fmpz.h"
#include "fixed.h"

/* Check fixed_rsqrt_ui_newton, fixed_rsqrt_newton and
   fixed_sqrt_newton against basecase references carrying XTRA extra
   limbs; the error bounds 4 B^-n / sqrt(a) translate, with
   a >= a_{an-1} / B, to 5 B^(XTRA+1) / sqrt(a_{an-1}). */

#undef XTRA
#undef check
#define XTRA 8
#define check newton_sqrt_check

static int
newton_sqrt_check(nn_srcptr got, slong gn, nn_srcptr ref, const fmpz_t bound)
{
    fmpz_t g, r;
    int ok;
    fmpz_init(g);
    fmpz_init(r);
    fmpz_set_ui_array(g, got, gn);
    fmpz_mul_2exp(g, g, FLINT_BITS * XTRA);
    fmpz_set_ui_array(r, ref, gn + XTRA);
    fmpz_sub(g, g, r);
    fmpz_abs(g, g);
    ok = (fmpz_cmp(g, bound) <= 0);
    fmpz_clear(g);
    fmpz_clear(r);
    return ok;
}

TEST_FUNCTION_START(fixed_sqrt_newton, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong n, an;
        nn_ptr A, Q, R;
        ulong a1;
        fmpz_t bound, t;

        if (n_randint(state, 8))
        {
            n = 1 + n_randint(state, 50);
            an = 1 + n_randint(state, 50);
        }
        else
        {
            n = 1 + n_randint(state, 600);
            an = 1 + n_randint(state, 600);
        }

        A = flint_malloc(an * sizeof(ulong));
        Q = flint_malloc((n + 2) * sizeof(ulong));
        R = flint_malloc((n + XTRA + 2) * sizeof(ulong));
        fmpz_init(bound);
        fmpz_init(t);

        flint_mpn_rrandom(A, state, an);
        if (n_randint(state, 5) == 0)
            A[an - 1] = 1 + n_randint(state, 100);
        if (A[an - 1] == 0)
            A[an - 1] = 1;

        fmpz_set_ui(t, A[an - 1]);
        fmpz_sqrt(t, t);
        fmpz_add_ui(t, t, 1);
        fmpz_set_ui(bound, 5);
        fmpz_mul_2exp(bound, bound, FLINT_BITS * (XTRA + 1));
        fmpz_tdiv_q(bound, bound, t);
        fmpz_add_ui(bound, bound, 8);

        fixed_rsqrt_newton(Q, A, an, n);
        fixed_rsqrt_newton_basecase(R, A, an, n + XTRA);

        if (!check(Q, n + 2, R, bound))
            TEST_FUNCTION_FAIL("rsqrt: n = %wd, an = %wd\n", n, an);

        fixed_sqrt_newton(Q, A, an, n);
        fixed_sqrt_newton_rsqrtmul(R, A, an, n + XTRA);

        fmpz_mul_ui(bound, bound, 2);
        if (!check(Q, n + 2, R, bound))
            TEST_FUNCTION_FAIL("sqrt: n = %wd, an = %wd\n", n, an);

        a1 = 2 + n_randint(state, UWORD_MAX - 2);
        fmpz_set_ui(bound, 3);
        fmpz_mul_2exp(bound, bound, FLINT_BITS * XTRA);

        fixed_rsqrt_ui_newton(Q, a1, n);
        fixed_rsqrt_ui_newton_basecase(R, a1, n + XTRA);

        if (!check(Q, n, R, bound))
            TEST_FUNCTION_FAIL("rsqrt_ui: n = %wd, a = %wu\n", n, a1);

        flint_free(A);
        flint_free(Q);
        flint_free(R);
        fmpz_clear(bound);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
