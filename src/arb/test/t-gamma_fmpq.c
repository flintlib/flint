/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq.h"
#include "arb.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

TEST_FUNCTION_START(arb_gamma_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t r, s;
        fmpq_t q;
        slong accuracy, prec, pp, qq;

        prec = 2 + n_randint(state, 1 << n_randint(state, 12));
        prec += 20;

        arb_init(r);
        arb_init(s);
        fmpq_init(q);

        pp = -100 + n_randint(state, 10000);
        qq = 1 + n_randint(state, 20);
        fmpq_set_si(q, pp, qq);

        arb_gamma_fmpq(r, q, prec);

        arb_set_fmpq(s, q, prec);
        arb_gamma(s, s, prec);

        if (!arb_overlaps(r, s))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("q = "); fmpq_print(q); flint_printf("\n\n");
            flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
            flint_printf("s = "); arb_printd(s, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        if (!(fmpz_is_one(fmpq_denref(q)) && fmpz_sgn(fmpq_numref(q)) <= 0)
            && FLINT_ABS(pp / qq) < 10)
        {
            accuracy = arb_rel_accuracy_bits(r);

            if (accuracy < prec - 6)
            {
                flint_printf("FAIL: poor accuracy\n\n");
                flint_printf("prec = %wd\n", prec);
                flint_printf("q = "); fmpq_print(q); flint_printf("\n\n");
                flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(r);
        arb_clear(s);
        fmpq_clear(q);
    }

    for (iter = 0; iter < 50 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t r, s;
        fmpq_t q;
        slong accuracy, prec;

        prec = 2 + n_randint(state, 25000);

        arb_init(r);
        arb_init(s);
        fmpq_init(q);

        fmpz_randtest(fmpq_numref(q), state, 3 + n_randlimb(state) % 30);
        fmpz_randtest_not_zero(fmpq_denref(q), state, 3 + n_randlimb(state) % 30);
        fmpq_canonicalise(q);

        arb_gamma_fmpq(r, q, prec);

        arb_set_fmpq(s, q, prec);
        arb_gamma(s, s, prec);

        if (!arb_overlaps(r, s))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("q = "); fmpq_print(q); flint_printf("\n\n");
            flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
            flint_printf("s = "); arb_printd(s, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        if (!(fmpz_is_one(fmpq_denref(q)) && fmpz_sgn(fmpq_numref(q)) <= 0) &&
            fabs(fmpq_get_d(q)) < 10.0)
        {
            accuracy = arb_rel_accuracy_bits(r);

            if (accuracy < prec - 6)
            {
                flint_printf("FAIL: poor accuracy\n\n");
                flint_printf("prec = %wd\n", prec);
                flint_printf("q = "); fmpq_print(q); flint_printf("\n\n");
                flint_printf("r = "); arb_printn(r, prec / 3.33, 0); flint_printf("\n\n");
                flint_printf("s = "); arb_printn(s, prec / 3.33, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(r);
        arb_clear(s);
        fmpq_clear(q);
    }

    TEST_FUNCTION_END(state);
}
