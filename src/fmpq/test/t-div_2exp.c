/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_div_2exp, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpq_t a, b;
        mpq_t c, d, e;
        flint_bitcnt_t exp;
        int aliasing;

        fmpq_init(a);
        fmpq_init(b);

        mpq_init(c);
        mpq_init(d);
        mpq_init(e);

        fmpq_randtest(a, state, 200);
        switch (n_randint(state, 5))
        {
            case 0:
                fmpz_mul_2exp(fmpq_numref(a), fmpq_numref(a), n_randint(state, 200));
                fmpq_canonicalise(a);
                break;
            case 1:
                fmpz_mul_2exp(fmpq_denref(a), fmpq_denref(a), n_randint(state, 200));
                fmpq_canonicalise(a);
                break;
        }
        fmpq_get_mpq(c, a);

        exp = n_randint(state, 200);

        aliasing = n_randint(state, 2);

        if (aliasing == 0)
        {
            fmpq_div_2exp(b, a, exp);
        }
        else
        {
            fmpq_set(b, a);
            fmpq_div_2exp(b, b, exp);
        }

        mpq_div_2exp(d, c, exp);

        fmpq_get_mpq(e, b);

        result = (mpq_cmp(d, e) == 0) && fmpq_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Qd, d = %Qd, e = %Qd, exp = %Md\n", c, d, e, exp);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);

        mpq_clear(c);
        mpq_clear(d);
        mpq_clear(e);
    }

    TEST_FUNCTION_END(state);
}
