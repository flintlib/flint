/*
    Copyright (C) 2010-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

void numerical_test(fmpq_t res, slong n, double ans)
{
    const double tol = 1e-13;
    double val, err;

    fmpq_harmonic_ui(res, n);
    val = fmpq_get_d(res);
    err = FLINT_ABS(val - ans);

    if (err > tol)
    {
        flint_printf("FAIL: %wd %.16f %.16f\n", n, val, ans);
        fflush(stdout);
        flint_abort();
    }
}

TEST_FUNCTION_START(fmpq_harmonic_ui, state)
{
    ulong i;
    fmpq_t s, t, u;

    fmpq_init(s);
    fmpq_init(t);
    fmpq_init(u);

    for (i = 0; i < 1000; i++)
    {
        if (i > 0)
        {
            fmpq_set_si(u, 1, i);
            fmpq_add(s, s, u);
        }

        fmpq_harmonic_ui(t, i);

        if (!fmpq_equal(t, s))
        {
            flint_printf("FAIL: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }
    }

    numerical_test(t, 1000, 7.4854708605503449127);
    numerical_test(t, 1001, 7.4864698615493459117);
    numerical_test(t, 1002, 7.4874678655413618797);
    numerical_test(t, 1003, 7.4884648745144426375);

    numerical_test(t, 10000, 9.7876060360443822642);
    numerical_test(t, 10001, 9.7877060260453821642);
    numerical_test(t, 10002, 9.7878060060493813643);
    numerical_test(t, 10003, 9.7879059760583786652);
    numerical_test(t, 10004, 9.7880059360743722677);

    numerical_test(t, 100000, 12.090146129863427947);

    if (flint_test_multiplier() > 10)
    {
        numerical_test(t, 20000, 10.480728217229327573);
        numerical_test(t, 30000, 10.886184992119899362);
        numerical_test(t, 40000, 11.173862897945522882);
        numerical_test(t, 50000, 11.397003949278482638);
        numerical_test(t, 60000, 11.579323839415955783);
        numerical_test(t, 70000, 11.733473328773164956);
        numerical_test(t, 80000, 11.867003828544530692);
        numerical_test(t, 90000, 11.984786169759202469);

        numerical_test(t, 100001, 12.090156129763428947);
        numerical_test(t, 100002, 12.090166129563432947);
        numerical_test(t, 100003, 12.090176129263441947);
        numerical_test(t, 100004, 12.090186128863457946);

        numerical_test(t, 300000, 13.188755085205611713);
        numerical_test(t, 500000, 13.699580042305528322);
        numerical_test(t, 700000, 14.036051993212618803);
        numerical_test(t, 900000, 14.287366262763433338);
    }

    fmpq_clear(s);
    fmpq_clear(t);
    fmpq_clear(u);

    TEST_FUNCTION_END(state);
}
