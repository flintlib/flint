/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_aors_ui, state)
{
    int i, result;

    for (i = 0; i < 20000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g, tst;
        mpz_t mf, mg;
        ulong x;
        int type;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(tst);

        mpz_init(mf);
        mpz_init(mg);

        type = n_randint(state, 5);
        if (type < 4)
        {
            fmpz_randtest(g, state, 200);
            x = n_randtest(state);
        }
        else
        {
            fmpz_randbits(g, state, FLINT_BITS + 1);
            x = n_randbits(state, FLINT_BITS);
        }

        fmpz_get_mpz(mg, g);

        switch (type)
        {
            /* Add */
            case 0:
                fmpz_add_ui(f, g, x);
                break;

            /* Add, aliased */
            case 1:
                fmpz_set(f, g);
                fmpz_add_ui(f, f, x);
                break;

            /* Sub */
            case 2:
                fmpz_sub_ui(f, g, x);
                break;

            /* Sub, aliased */
            case 3:
                fmpz_set(f, g);
                fmpz_sub_ui(f, f, x);
                break;

            /* size(g) = 2 but f may be small */
            case 4:
                if (fmpz_sgn(g) < 0)
                    fmpz_add_ui(f, g, x);
                else
                    fmpz_sub_ui(f, g, x);
                break;

            default: FLINT_UNREACHABLE;
        }

        if (type < 2 || (type == 4 && fmpz_sgn(g) < 0))
            flint_mpz_add_ui(mf, mg, x);
        else
            flint_mpz_sub_ui(mf, mg, x);

        fmpz_set_mpz(tst, mf);

        result = fmpz_equal(f, tst) && _fmpz_is_canonical(f);

        if (!result)
            flint_throw(FLINT_TEST_FAIL,
                    "f = %{fmpz}\n"
                    "g = %{fmpz}\n"
                    "x = %wu\n\n"
                    "Correct result via GMP: %{fmpz}\n",
                    f, g, x, tst);

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(tst);

        mpz_clear(mf);
        mpz_clear(mg);
    }

    TEST_FUNCTION_END(state);
}
