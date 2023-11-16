/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpzi.h"

TEST_FUNCTION_START(fmpzi_divrem_approx, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpzi_t x, y, q, r, s, t;
        fmpz_t ny, nr;
        int aliasing;

        fmpzi_init(x);
        fmpzi_init(y);
        fmpzi_init(q);
        fmpzi_init(r);
        fmpzi_init(s);
        fmpzi_init(t);
        fmpz_init(ny);
        fmpz_init(nr);

        fmpzi_randtest(x, state, 2 + n_randint(state, 600));
        do {
            fmpzi_randtest(y, state, 2 + n_randint(state, 600));
        } while (fmpzi_is_zero(y));

        fmpzi_randtest(q, state, 2 + n_randint(state, 600));
        fmpzi_randtest(r, state, 2 + n_randint(state, 600));

        if (n_randint(state, 3) == 0)
        {
            fmpzi_mul(x, y, q);
            if (n_randint(state, 3) == 0)
                fmpzi_add(x, x, r);
        }

        aliasing = n_randint(state, 7);
        switch (aliasing)
        {
            case 0:
                fmpzi_divrem_approx(q, r, x, y);
                break;
            case 1:
                fmpzi_set(q, x);
                fmpzi_divrem_approx(q, r, q, y);
                break;
            case 2:
                fmpzi_set(q, y);
                fmpzi_divrem_approx(q, r, x, q);
                break;
            case 3:
                fmpzi_set(r, x);
                fmpzi_divrem_approx(q, r, r, y);
                break;
            case 4:
                fmpzi_set(r, y);
                fmpzi_divrem_approx(q, r, x, r);
                break;
            case 5:
                fmpzi_set(q, x);
                fmpzi_set(r, y);
                fmpzi_divrem_approx(q, r, q, r);
                break;
            default:
                fmpzi_set(q, y);
                fmpzi_set(r, x);
                fmpzi_divrem_approx(q, r, r, q);
                break;
        }

        fmpzi_mul(s, q, y);
        fmpzi_add(t, s, r);

        fmpzi_norm(ny, y);
        fmpzi_norm(nr, r);

        if (!fmpzi_equal(t, x) || !(fmpz_cmp(nr, ny) < 0))
        {
            flint_printf("FAIL\n");
            flint_printf("aliasing = %d\n", aliasing);
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("q = "); fmpzi_print(q); printf("\n");
            flint_printf("r = "); fmpzi_print(r); printf("\n");
            flint_printf("s = "); fmpzi_print(s); printf("\n");
            flint_printf("t = "); fmpzi_print(t); printf("\n");
            flint_printf("nr = "); fmpz_print(nr); printf("\n");
            flint_printf("ny = "); fmpz_print(ny); printf("\n");
            flint_abort();
        }

        fmpzi_clear(x);
        fmpzi_clear(y);
        fmpzi_clear(q);
        fmpzi_clear(r);
        fmpzi_clear(s);
        fmpzi_clear(t);
        fmpz_clear(ny);
        fmpz_clear(nr);
    }

    TEST_FUNCTION_END(state);
}
