/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_ui_vec_prod, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nn_ptr a = NULL, b = NULL;
        fmpz_t x, y, z;
        slong i;

        slong len1 = n_randint(state, 100);
        slong len2 = n_randint(state, 100);

        a = flint_realloc(a, sizeof(ulong) * (len1 + len2));
        b = a + len1;

        for (i = 0; i < len1 + len2; i++)
            a[i] = n_randtest(state);

        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);

        _fmpz_ui_vec_prod(x, a, len1);
        _fmpz_ui_vec_prod(y, b, len2);
        fmpz_mul(x, x, y);
        _fmpz_ui_vec_prod(z, a, len1 + len2);

        result = (fmpz_equal(x, z));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("%ld ", len1);
            for (i = 0; i < len1; i++)
                flint_printf("%lu ", a[i]);
            flint_printf("\n\n");
            flint_printf("%ld ", len2);
            for (i = 0; i < len2; i++)
                flint_printf("%lu ", b[i]);
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(a);

        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
    }

    TEST_FUNCTION_END(state);
}
