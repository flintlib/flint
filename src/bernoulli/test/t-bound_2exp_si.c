/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arf.h"
#include "bernoulli.h"

double log2bern_approx(double n)
{
    return 1 + ((n+0.5)*log(n) - n - (n-0.5)*log(2*3.14159265358979323)) * (1. / log(2));
}

TEST_FUNCTION_START(bernoulli_bound_2exp_si, state)
{
    slong i, bound;
    double a, b;
    fmpq_t q;
    arf_t t;

    fmpq_init(q);
    arf_init(t);

    for (i = 0; i < 1000; i++)
    {
        arith_bernoulli_number(q, i);
        bound = bernoulli_bound_2exp_si(i);

        arf_set_round_fmpz(t, fmpq_numref(q), 32, ARF_RND_UP);
        arf_div_fmpz(t, t, fmpq_denref(q), 32, ARF_RND_UP);

        if (arf_cmpabs_2exp_si(t, bound) > 0)
        {
            flint_printf("FAIL: %wd\n", i);
            arf_print(t); flint_printf("\n\n");
            flint_printf("%wd\n", bound); flint_printf("\n\n");
            flint_abort();
        }
    }

    fmpq_clear(q);
    arf_clear(t);

    for (i = 100; i < 4000000; i += 1)
    {
        i += (i & 1);
        a = bernoulli_bound_2exp_si(i);
        b = log2bern_approx(i);

        if (a < b || a > 1.01 * b)
        {
            flint_printf("FAIL: %wd\n", i);
            flint_printf("%wd: %f %f %f\n", i, a, b, (float) a / b);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
