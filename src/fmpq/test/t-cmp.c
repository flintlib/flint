/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "fmpz.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_cmp, state)
{
    int i;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpq_t x, y;
        mpq_t X, Y;
        int c1, c2;
        int type;

        fmpq_init(x);
        fmpq_init(y);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);

        type = n_randint(state, 4);

        switch (type)
        {
            case 0:
                fmpq_randtest(y, state, 200);
                c1 = fmpq_cmp(x, y);
                break;
            case 1:
                fmpz_randtest(fmpq_denref(y), state, 200);
                c1 = fmpq_cmp_fmpz(x, fmpq_numref(y));
                break;
            case 2:
            {
                slong ys = z_randtest(state);
                c1 = fmpq_cmp_si(x, ys);
                fmpz_set_si(fmpq_numref(y), ys);
                break;
            }
            case 3:
            {
                ulong ys = n_randtest(state);
                c1 = fmpq_cmp_ui(x, ys);
                fmpz_set_ui(fmpq_numref(y), ys);
                break;
            }
            default: FLINT_UNREACHABLE;
        }

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);
        c2 = mpq_cmp(X, Y);

        if (FLINT_SGN(c1) != FLINT_SGN(c2))
        {
            flint_throw(FLINT_TEST_FAIL,
                    "x = %{fmpq}\n"
                    "y = %{fmpq}\n"
                    "%s(x, y) = %d\n"
                    "cmp(X, Y) = %d\n",
                    x, y,
                    type == 0 ? "cmp" : type == 1 ? "cmp_fmpz" : type == 2 ? "cmp_si" : "cmp_ui",
                    c1,
                    c2);
        }

        fmpq_clear(x);
        fmpq_clear(y);

        mpq_clear(X);
        mpq_clear(Y);
    }

    TEST_FUNCTION_END(state);
}
