/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_frexp, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z;
        fmpz_t e, f;
        int ok;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        fmpz_init(e);
        fmpz_init(f);

        arf_randtest_special(x, state, 200, 100);
        arf_randtest_special(y, state, 200, 100);
        fmpz_randtest(e, state, 100);
        fmpz_randtest(f, state, 100);

        arf_frexp(y, e, x);
        arf_mul_2exp_fmpz(z, y, e);

        ok = 1;

        if (!arf_equal(z, x)) ok = 0;
        if (arf_is_special(x) && !fmpz_is_zero(e)) ok = 0;
        if (!arf_is_special(x) && !(arf_cmpabs_2exp_si(y, 0) < 0
            && arf_cmpabs_2exp_si(y, -1) >= 0)) ok = 0;

        if (!ok)
        {
            printf("FAIL\n");
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = "); arf_print(y); printf("\n\n");
            printf("z = "); arf_print(z); printf("\n\n");
            printf("e = "); fmpz_print(e); printf("\n\n");
            printf("f = "); fmpz_print(f); printf("\n\n");
            flint_abort();
        }

        arf_frexp(x, f, x);

        if (!arf_equal(x, y) || !fmpz_equal(e, f))
        {
            printf("FAIL (aliasing)\n");
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = "); arf_print(y); printf("\n\n");
            printf("z = "); arf_print(z); printf("\n\n");
            printf("e = "); fmpz_print(e); printf("\n\n");
            printf("f = "); fmpz_print(f); printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        fmpz_clear(e);
        fmpz_clear(f);
    }

    TEST_FUNCTION_END(state);
}
