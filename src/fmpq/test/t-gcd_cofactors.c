/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_gcd_cofactors, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_t g, a, b, g1, t1, t2;
        fmpz_t abar, bbar, abar1, bbar1;

        fmpq_init(g);
        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(g1);
        fmpq_init(t1);
        fmpq_init(t2);
        fmpz_init(abar);
        fmpz_init(bbar);
        fmpz_init(abar1);
        fmpz_init(bbar1);

        fmpq_randtest(g, state, 200);
        fmpq_randtest(a, state, 200);
        fmpq_randtest(b, state, 200);

        fmpq_gcd_cofactors(g, abar, bbar, a, b);
        fmpq_gcd(g1, a, b);

        if (!fmpq_equal(g, g1))
        {
            flint_printf("FAIL: check gcd\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_mul_fmpz(t1, g, abar);
        fmpq_mul_fmpz(t2, g, bbar);
        if (!fmpq_equal(t1, a) || !fmpq_equal(t2, b))
        {
            flint_printf("FAIL: check cofactors\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_set(t1, a);
        fmpq_gcd_cofactors(t1, abar1, bbar1, t1, b);
        if (!fmpq_equal(t1, g) ||
            !fmpz_equal(abar1, abar) ||!fmpz_equal(bbar1, bbar))
        {
            flint_printf("FAIL: check aliasing first argument\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_set(t1, b);
        fmpq_gcd_cofactors(t1, abar1, bbar1, a, t1);
        if (!fmpq_equal(t1, g) ||
            !fmpz_equal(abar1, abar) ||!fmpz_equal(bbar1, bbar))
        {
            flint_printf("FAIL: check aliasing second argument\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(g);
        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(g1);
        fmpq_clear(t1);
        fmpq_clear(t2);
        fmpz_clear(abar);
        fmpz_clear(bbar);
        fmpz_clear(abar1);
        fmpz_clear(bbar1);
    }

    TEST_FUNCTION_END(state);
}
