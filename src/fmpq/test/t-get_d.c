/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_get_d, state)
{
    const flint_bitcnt_t bound = 1000;
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        flint_bitcnt_t num_bits, den_bits;
        int sgn, result;
        double d;
        fmpq_t a;

        fmpq_init(a);

        fmpq_randtest(a, state, 9999);

        num_bits = fmpz_bits(fmpq_numref(a));
        den_bits = fmpz_bits(fmpq_denref(a));

        if (   (num_bits >= den_bits && num_bits - den_bits < bound)
            || (den_bits >= num_bits && den_bits - num_bits < bound))
        {
            /* the exponent range is such that d should be finite */
            d = fmpq_get_d(a);
            sgn = fmpq_sgn(a);
            result = (sgn == 0) ? (d == 0) : (sgn < 0) ? (d < 0) : (d > 0);

            if (!result)
            {
                flint_printf("FAIL:\ncheck sign of result matches\n");
                flint_printf("a = "); fmpq_print(a); flint_printf("\n");
                flint_printf("num_bits = %wu\n", num_bits);
                flint_printf("den_bits = %wu\n", den_bits);
                flint_printf("d = %f\n", d);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_clear(a);
    }

    TEST_FUNCTION_END(state);
}
