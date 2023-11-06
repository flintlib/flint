/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_add_fmpz, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpq_t a;
        fmpz_t b;
        fmpq_t c, d, e, f;
        int aliasing;

        fmpq_init(a);
        fmpz_init(b);
        fmpq_init(c);
        fmpq_init(d);
        fmpq_init(e);
        fmpq_init(f);

        fmpq_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        fmpq_set(d, a);
        fmpq_set_fmpz(e, b);

        aliasing = n_randint(state, 4);

        if (aliasing == 0)
        {
            fmpq_add_fmpz(c, a, b);
        }
        else if (aliasing == 1)
        {
            fmpq_set_fmpz(e, fmpq_numref(a));
            fmpq_add_fmpz(c, a, fmpq_numref(a));
        }
        else if (aliasing == 2)
        {
            fmpq_set(c, a);
            fmpq_add_fmpz(c, c, b);
        }
        else
        {
            fmpq_set_fmpz(c, b);
            fmpq_add_fmpz(c, a, fmpq_numref(c));
        }

        fmpq_add(f, d, e);

        result = (fmpq_cmp(f, c) == 0) && fmpq_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            printf("c = "); fmpq_print(c);
            printf(", d = "); fmpq_print(d);
            printf(", e = "); fmpq_print(e);
            printf(", f = "); fmpq_print(f); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpz_clear(b);
        fmpq_clear(c);
        fmpq_clear(d);
        fmpq_clear(e);
        fmpq_clear(f);
    }

    TEST_FUNCTION_END(state);
}
