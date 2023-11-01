/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_get_nmod_poly, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t A;
        nmod_poly_t M, M2;
        slong length;
        mp_limb_t mod;

        length = n_randint(state, 50);

        mod = n_randtest_prime(state, 0);

        nmod_poly_init(M, mod);
        nmod_poly_init(M2, mod);
        fmpz_poly_init(A);

        nmod_poly_randtest(M, state, length);

        if (i % 2 == 0)
            fmpz_poly_set_nmod_poly(A, M);
        else
            fmpz_poly_set_nmod_poly_unsigned(A, M);

        fmpz_poly_scalar_mul_ui(A, A, UWORD(2));
        nmod_poly_add(M, M, M);
        fmpz_poly_get_nmod_poly(M2, A);

        if (!nmod_poly_equal(M, M2))
        {
            flint_printf("FAIL!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(A);
        nmod_poly_clear(M);
        nmod_poly_clear(M2);
    }

    TEST_FUNCTION_END(state);
}
