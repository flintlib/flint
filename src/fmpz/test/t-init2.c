/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_init2, state)
{
    int i;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;

        fmpz_init2(a, n_randint(state, 100));
        fmpz_clear(a);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;

        fmpz_init(a);
        fmpz_randtest(a, state, SMALL_FMPZ_BITCOUNT_MAX);

        _fmpz_promote_val(a);
        _fmpz_demote_val(a);

        fmpz_clear(a);
    }

    TEST_FUNCTION_END(state);
}
