/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly_q.h"

TEST_FUNCTION_START(fmpz_poly_q_init_clear, state)
{
    int i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a;
        slong len1 = n_randint(state, 50);
        slong len2 = n_randint(state, 50);
        flint_bitcnt_t bits1 = n_randint(state, 50);
        flint_bitcnt_t bits2 = n_randint(state, 50);

        fmpz_poly_q_init(a);
        fmpz_poly_q_randtest(a, state, len1, bits1, len2, bits2);
        fmpz_poly_q_clear(a);
    }

    TEST_FUNCTION_END(state);
}
