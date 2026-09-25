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
#include "ulong_extras.h"
#include "arb.h"
#include "fixed.h"

/* fixed_neglog_newton and fixed_atan_newton agree with arb within
   their budgets over the whole input ranges (random inputs, inputs
   with long runs of zero or one bits including the endpoints, and
   inputs whose top limbs are zero so that the starting value is
   degenerate), for the tunable parameters,
   including sizes that exercise the recursive starting value and
   the windowed products at very small sizes */

static double
_newton_ulps(nn_srcptr y, slong n, const arb_t exact)
{
    arb_t va, d;
    fmpz_t f;
    mag_t mb;
    double u;

    arb_init(va); arb_init(d); fmpz_init(f); mag_init(mb);
    fmpz_set_ui_array(f, y, n);
    arb_set_fmpz(va, f);
    arb_mul_2exp_si(va, va, -FLINT_BITS * n);
    arb_sub(d, va, exact, ARF_PREC_EXACT);
    arb_get_mag(mb, d);
    mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
    u = mag_get_d(mb);
    arb_clear(va); arb_clear(d); fmpz_clear(f); mag_clear(mb);
    return u;
}

/* random n-limb fraction with structured bit patterns some of the time */
static void
_randlimbs(nn_ptr x, slong n, flint_rand_t state)
{
    slong k;

    switch (n_randint(state, 6))
    {
        case 0:
            flint_mpn_urandomb(x, state, FLINT_BITS * n);
            for (k = 0; k < n; k++)
                if (n_randint(state, 3) == 0)
                    x[k] = n_randint(state, 2) ? ~UWORD(0) : 0;
            break;
        case 1:
            flint_mpn_zero(x, n);
            k = n_randint(state, n);
            flint_mpn_urandomb(x, state, FLINT_BITS * (k + 1));
            break;
        case 2:
            for (k = 0; k < n; k++)
                x[k] = ~UWORD(0);
            k = n_randint(state, n);
            x[k] = n_randtest(state);
            break;
        case 3:
            flint_mpn_zero(x, n);
            x[n_randint(state, n)] = UWORD(1) << n_randint(state, FLINT_BITS);
            break;
        default:
            flint_mpn_urandomb(x, state, FLINT_BITS * n);
    }
}

TEST_FUNCTION_START(fixed_newton, state)
{
    slong iter;

    for (iter = 0; iter < 200 + 100 * flint_test_multiplier(); iter++)
    {
        slong n, N;
        int forward, which;
        nn_ptr x, y;
        arb_t xa, r;
        fmpz_t f;
        double u, bound;

        if (iter % 19 == 0)
            n = 1 + n_randint(state, 3000);   /* through the recursive start */
        else if (iter % 5 == 0)
            n = 1 + n_randint(state, 300);
        else
            n = 1 + n_randint(state, 40);

        which = n_randint(state, 2);
        forward = (n <= 200) ? n_randint(state, 3) : 2 * n_randint(state, 2);
        N = n_randint(state, 3) ? 0 : 1 + n_randint(state, 16);

        x = flint_malloc(n * sizeof(ulong));
        y = flint_malloc(n * sizeof(ulong));
        arb_init(xa); arb_init(r); fmpz_init(f);

        _randlimbs(x, n, state);
        if (which == 0)
            x[n - 1] |= UWORD(1) << (FLINT_BITS - 1);   /* [1/2, 1) */

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);

        if (which == 0)
        {
            arb_log(r, xa, FLINT_BITS * n + 64);
            arb_neg(r, r);
            _fixed_neglog_newton_tune(y, x, n, forward, N);
            bound = FIXED_NEGLOG_NEWTON_MAX_ERR;
        }
        else
        {
            arb_atan(r, xa, FLINT_BITS * n + 64);
            _fixed_atan_newton_tune(y, x, n, forward, N);
            bound = FIXED_ATAN_NEWTON_MAX_ERR;
        }

        u = _newton_ulps(y, n, r);
        if (!(u <= bound))
            TEST_FUNCTION_FAIL("%s: n = %wd, forward = %d, N = %wd: "
                "%g ulps\nx = %{ulong*}\ny = %{ulong*}\n",
                which ? "atan" : "neglog", n, forward, N, u, x, n, y, n);

        /* the public entries */
        if (which == 0)
            fixed_neglog_newton(y, x, n);
        else
            fixed_atan_newton(y, x, n);
        u = _newton_ulps(y, n, r);
        if (!(u <= bound))
            TEST_FUNCTION_FAIL("%s (default): n = %wd: %g ulps\nx = %{ulong*}\n",
                which ? "atan" : "neglog", n, u, x, n);

        /* the AGM logarithm: the same budget, and its ball contains
           the truth */
        if (which == 0)
        {
            fball_t rb;
            arb_t ra;
            fball_init(rb);
            arb_init(ra);
            _fixed_neglog_agm_tune(y, x, n, N);
            u = _newton_ulps(y, n, r);
            if (!(u <= bound))
                TEST_FUNCTION_FAIL("neglog_agm: n = %wd, N = %wd: %g ulps\nx = %{ulong*}\n",
                    n, N, u, x, n);
            _fball_neglog_agm(rb, x, n, N);
            fball_get_arb(ra, rb);
            if (!arb_overlaps(ra, r))
                TEST_FUNCTION_FAIL("neglog_agm (ball): n = %wd, N = %wd\nx = %{ulong*}\n",
                    n, N, x, n);
            fball_clear(rb);
            arb_clear(ra);
        }

        /* the ball worker contains the truth */
        {
            fball_t rb;
            arb_t ra;
            fball_init(rb);
            arb_init(ra);
            if (which == 0)
                _fball_neglog_newton(rb, x, n, forward, N);
            else
                _fball_atan_newton(rb, x, n, forward, N);
            fball_get_arb(ra, rb);
            if (!arb_overlaps(ra, r))
                TEST_FUNCTION_FAIL("%s (ball): n = %wd, forward = %d, N = %wd: "
                    "ball does not contain the value\nx = %{ulong*}\n",
                    which ? "atan" : "neglog", n, forward, N, x, n);
            fball_clear(rb);
            arb_clear(ra);
        }

        flint_free(x); flint_free(y);
        arb_clear(xa); arb_clear(r); fmpz_clear(f);
    }

    TEST_FUNCTION_END(state);
}
