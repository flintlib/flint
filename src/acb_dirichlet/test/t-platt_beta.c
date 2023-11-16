/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

static void
_arb_div_ui_ui(arb_t res, ulong a, ulong b, slong prec)
{
    arb_set_ui(res, a);
    arb_div_ui(res, res, b, prec);
}

static int
_arb_lt_d(const arb_t a, double d)
{
    int result;
    arb_t x;
    arb_init(x);
    arb_set_d(x, d);
    result = arb_lt(a, x);
    arb_clear(x);
    return result;
}

TEST_FUNCTION_START(acb_dirichlet_platt_beta, state)
{
    slong iter;
    arb_t x, t, t0, expe;

    arb_init(x);
    arb_init(t);
    arb_init(t0);
    arb_init(expe);

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t beta, a, b, c;
        acb_t z;
        slong prec, mbits;

        prec = 2 + n_randint(state, 300);
        mbits = 2 + n_randint(state, 20);
        arb_randtest(t, state, prec, mbits);
        arb_randtest(t0, state, prec, mbits);
        arb_abs(t, t);
        arb_abs(t0, t0);
        arb_add(x, t, t0, prec);
        arb_const_e(expe, prec);
        arb_exp(expe, expe, prec);

        if (!arb_is_nonnegative(t) || !arb_gt(t0, expe) || !_arb_lt_d(x, 1e8))
        {
            continue;
        }

        arb_init(beta);
        arb_init(a);
        arb_init(b);
        arb_init(c);
        acb_init(z);

        acb_dirichlet_platt_beta(beta, t0, prec);

        /* Lemma A.10 in "Isolating some non-trivial zeros of zeta" */
        arb_pow(a, x, beta, prec);
        arb_mul_ui(a, a, 3, prec);
        acb_dirichlet_platt_scaled_lambda(c, t, prec);
        arb_abs(c, c);
        if (arb_gt(c, a))
        {
            flint_printf("FAIL: Lemma A.10 |f(t)|\n\n");
            flint_printf("iter = %wd\n\n", iter);
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("t = "); arb_printd(t, 15); flint_printf("\n\n");
            flint_printf("t0 = "); arb_printd(t0, 15); flint_printf("\n\n");
            flint_abort();
        }

        /* In Lemma A.10 proof in "Isolating some non-trivial zeros of zeta" */
        arb_pow(a, x, beta, prec);
        _arb_div_ui_ui(b, 732, 1000, prec);
        arb_mul(a, a, b, prec);
        acb_set_d(z, 0.5);
        arb_set(acb_imagref(z), x);
        acb_zeta(z, z, prec);
        acb_abs(c, z, prec);
        if (arb_gt(c, a))
        {
            flint_printf("FAIL: Lemma A.10 |zeta(1/2 + i(t + t0))|\n\n");
            flint_printf("iter = %wd\n\n", iter);
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("t = "); arb_printd(t, 15); flint_printf("\n\n");
            flint_printf("t0 = "); arb_printd(t0, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(beta);
        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        acb_clear(z);
    }

    arb_clear(x);
    arb_clear(t);
    arb_clear(t0);
    arb_clear(expe);

    TEST_FUNCTION_END(state);
}
