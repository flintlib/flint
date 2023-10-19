/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_init_set_readonly, state)
{
    int i;

    /* Create some small fmpq rationals, clear the mpq_t */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, SMALL_FMPZ_BITCOUNT_MAX);

        mpq_init(z);
        fmpq_get_mpq(z, f);

        {
            fmpq_t g;

            fmpq_init_set_readonly(g, z);
            fmpq_clear_readonly(g);
        }

        mpq_clear(z);
    }

    /* Create some small fmpq ratioals, do *not* clear the mpq_t */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, SMALL_FMPZ_BITCOUNT_MAX);

        mpq_init(z);
        fmpq_get_mpq(z, f);

        {
            fmpq_t g;

            fmpq_init_set_readonly(g, z);
        }

        mpq_clear(z);
    }

    /* Create some more fmpq rationals */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, 2 * FLINT_BITS);
        mpq_init(z);
        fmpq_get_mpq(z, f);

        {
            fmpq_t g, h;

            fmpq_init_set_readonly(g, z);
            fmpq_init(h);
            fmpq_set_mpq(h, z);

            if (!fmpq_equal(g, h))
            {
                flint_printf("FAIL:\n\n");
                flint_printf("g = "), fmpq_print(g), flint_printf("\n");
                flint_printf("h = "), fmpq_print(h), flint_printf("\n");
                gmp_printf("z = %Qd\n", z);
            }

            fmpq_clear_readonly(g);
            fmpq_clear(h);
        }

        fmpq_clear(f);
        mpq_clear(z);
    }

    TEST_FUNCTION_END(state);
}
