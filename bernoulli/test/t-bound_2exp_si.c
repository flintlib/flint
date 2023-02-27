/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "bernoulli.h"

double log2bern_approx(double n)
{
    return 1 + ((n+0.5)*log(n) - n - (n-0.5)*log(2*3.14159265358979323)) * (1. / log(2));
}

int main()
{
    slong i, bound;
    double a, b;
    fmpq_t q;
    fmpr_t t;

    flint_printf("bound_2exp_si....");
    fflush(stdout);

    fmpq_init(q);
    fmpr_init(t);

    for (i = 0; i < 1000; i++)
    {
        arith_bernoulli_number(q, i);
        bound = bernoulli_bound_2exp_si(i);

        fmpr_set_round_fmpz(t, fmpq_numref(q), 32, FMPR_RND_UP);
        fmpr_div_fmpz(t, t, fmpq_denref(q), 32, FMPR_RND_UP);

        if (fmpr_cmpabs_2exp_si(t, bound) > 0)
        {
            flint_printf("FAIL: %wd\n", i);
            fmpr_print(t); flint_printf("\n\n");
            flint_printf("%wd\n", bound); flint_printf("\n\n");
            flint_abort();
        }
    }

    fmpq_clear(q);
    fmpr_clear(t);

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

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

