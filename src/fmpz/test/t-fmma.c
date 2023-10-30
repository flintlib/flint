/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_fmma, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d, f, g;
        fmpz *aptr, *bptr, *cptr, *dptr;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(f);
        fmpz_init(g);

        fmpz_randtest(f, state, 200);
        fmpz_randtest(g, state, 200);

        if (n_randint(state, 2)) { fmpz_randtest(a, state, 200); aptr = a; }
        else { fmpz_set(a, f); aptr = f; }
        if (n_randint(state, 2)) { fmpz_randtest(b, state, 200); bptr = b; }
        else { fmpz_set(b, f); bptr = f; }
        if (n_randint(state, 2)) { fmpz_randtest(c, state, 200); cptr = c; }
        else { fmpz_set(c, f); cptr = f; }
        if (n_randint(state, 2)) { fmpz_randtest(d, state, 200); dptr = d; }
        else { fmpz_set(d, f); dptr = f; }

        fmpz_fmma(f, aptr, bptr, cptr, dptr);

        fmpz_mul(g, a, b);
        fmpz_addmul(g, c, d);

        if (!fmpz_equal(f, g) || !_fmpz_is_canonical(f))
        {
            flint_printf("FAIL:\n");
            fmpz_print(a); flint_printf("\n");
            fmpz_print(b); flint_printf("\n");
            fmpz_print(c); flint_printf("\n");
            fmpz_print(d); flint_printf("\n");
            fmpz_print(f); flint_printf("\n");
            fmpz_print(g); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(f);
        fmpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}
