/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "arith.h"

TEST_FUNCTION_START(arith_bernoulli_number_denom, state)
{
    fmpz_t s, t;
    slong n;


    fmpz_init(s);
    fmpz_init(t);

    for (n = 0; n < 1000; n++)
    {
        arith_bernoulli_number_denom(t, n);
        fmpz_addmul_ui(s, t, n_nth_prime(n+1));
    }

    fmpz_set_str(t, "34549631155954474103407159", 10);

    if (!fmpz_equal(s, t))
    {
        flint_printf("FAIL: Hash disagrees with known value\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_clear(s);
    fmpz_clear(t);

    TEST_FUNCTION_END(state);
}
