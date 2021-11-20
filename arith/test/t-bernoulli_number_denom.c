/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "fmpz.h"
#include "ulong_extras.h"

int main()
{
    fmpz_t s, t;
    slong n;

    FLINT_TEST_INIT(state);

    flint_printf("bernoulli_number_denom....");
    fflush(stdout);

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

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
