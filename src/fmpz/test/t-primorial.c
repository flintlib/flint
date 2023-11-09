/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_primorial, state)
{
    ulong k;
    fmpz_t x;
    fmpz_t y;

    fmpz_init(x);
    fmpz_init(y);
    fmpz_set_ui(y, 1);

    for (k = 0; k < 2000; k++)
    {
       fmpz_primorial(x, k);
       if (n_is_prime(k))
          fmpz_mul_ui(y, y, k);
       if (!fmpz_equal(x, y) || !_fmpz_is_canonical(x))
       {
          flint_printf("FAIL:\n");
          flint_printf("primorial of %wu disagrees with direct product\n", k);
          fmpz_print(x);
          flint_printf("\n");
          fflush(stdout);
          flint_abort();
       }
    }

    fmpz_clear(x);
    fmpz_clear(y);

    TEST_FUNCTION_END(state);
}
