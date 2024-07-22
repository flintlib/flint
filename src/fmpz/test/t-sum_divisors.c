/*
    Copyright (C) 2024 Matthias Gessinger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

void _fmpz_sum_divisors_naive(fmpz_t f, fmpz_t g)
{
    fmpz_t sqrt, temp1, temp2, temp3, out;

    fmpz_init(sqrt);
    fmpz_init(temp1);
    fmpz_init(temp2);
    fmpz_init(temp3);
    fmpz_init(out);

    int is_neg = fmpz_cmp_ui(g, 0) < 1;

    fmpz_abs(f, g);
    fmpz_sqrt(sqrt, f);

    fmpz_set_si(temp1, 1);
    while (fmpz_cmp(temp1, sqrt) != 1)
    {
        fmpz_tdiv_qr(temp3, temp2, f, temp1);

        if (fmpz_is_zero(temp2))
        {
            fmpz_add(out, out, temp3);

            if (fmpz_cmp(temp1, temp3) != 0)
            {
                fmpz_add(out, out, temp1);
            }
        }

        fmpz_add_si(temp1, temp1, 1);
    }

    if (is_neg)
    {
        fmpz_neg(f, out);
    }
    else
    {
        fmpz_set(f, out);
    }

    fmpz_clear(out);
    fmpz_clear(temp1);
    fmpz_clear(temp2);
    fmpz_clear(temp3);
    fmpz_clear(sqrt);
}

TEST_FUNCTION_START(fmpz_sum_divisors, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_randtest(a, state, 32);

        aliasing = n_randint(state, 2);

        /* The reference result */
        _fmpz_sum_divisors_naive(c, a);

        /* general function */
        if (aliasing == 1)
        {
            fmpz_sum_divisors(a, a);
            fmpz_set(b, a);
        }
        else
        {
            fmpz_sum_divisors(b, a);
        }

        result = (fmpz_cmp(c, b) == 0);

        if (!result)
        {
            TEST_FUNCTION_FAIL(
                "Expected: %{fmpz}\n"
                "Got: %{fmpz}\n",
                c, b);
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
