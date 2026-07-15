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

/* Check fixed_inv_newton and fixed_div_newton against basecase
   references carrying XTRA extra limbs: with the documented error
   bound 4 B^-n / a and the reference's own few ulps at B^-(n+XTRA),
   the difference scaled by B^XTRA must stay below 5 B^(XTRA+1)
   divided by the denominator's top limb (1/a < B / a_{an-1}). */

#define XTRA 8
#define check newton_div_check

static int
newton_div_check(nn_srcptr got, slong gn, nn_srcptr ref, const fmpz_t bound)
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

TEST_FUNCTION_START(fixed_div_newton, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong n, an, bn;
        nn_ptr A, Q, R;
        fmpz_t bound;

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
        bn = 1 + n_randint(state, an);

        A = flint_malloc(an * sizeof(ulong));
        Q = flint_malloc((n + 2) * sizeof(ulong));
        R = flint_malloc((n + XTRA + 2) * sizeof(ulong));
        fmpz_init(bound);

        flint_mpn_rrandom(A, state, an);
        if (n_randint(state, 5) == 0)
            A[an - 1] = 1 + n_randint(state, 100);
        if (A[an - 1] == 0)
            A[an - 1] = 1;

        fmpz_set_ui(bound, 5);
        fmpz_mul_2exp(bound, bound, FLINT_BITS * (XTRA + 1));
        fmpz_tdiv_q_ui(bound, bound, A[an - 1]);
        fmpz_add_ui(bound, bound, 8);

        fixed_inv_newton(Q, A, an, n);
        fixed_inv_newton_basecase(R, A, an, n + XTRA);

        if (!check(Q, n + 2, R, bound))
            TEST_FUNCTION_FAIL("inv: n = %wd, an = %wd\n", n, an);

        fixed_div_newton(Q, A, bn, A, an, n);
        fixed_div_newton_invmul(R, A, bn, A, an, n + XTRA);

        /* both sides approximate: double the allowance */
        fmpz_mul_ui(bound, bound, 2);
        if (!check(Q, n + 2, R, bound))
            TEST_FUNCTION_FAIL("div: n = %wd, an = %wd, bn = %wd\n",
                n, an, bn);

        flint_free(A);
        flint_free(Q);
        flint_free(R);
        fmpz_clear(bound);
    }

    TEST_FUNCTION_END(state);
}
