/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_fib, state)
{
    /* trivial test because functions are based on gr */
    {
        fmpz_t n;
        arb_t t;
        fmpz_init_set_ui(n, 10);
        arb_init(t);
        arb_fib_fmpz(t, n, 32);
        if (!arb_equal_si(t, 55))
            flint_abort();
        arb_fib_ui(t, 11, 32);
        if (!arb_equal_si(t, 89))
            flint_abort();
        arb_clear(t);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
