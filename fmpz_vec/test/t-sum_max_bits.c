/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("sum_max_bits....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz * a;
        fmpz_t max, sum, t;
        slong sum_bits, max_bits;
        slong len = n_randint(state, 300);

        a = _fmpz_vec_init(len);
        fmpz_init(max);
        fmpz_init(sum);
        fmpz_init(t);

        _fmpz_vec_randtest(a, state, len, 300);

        _fmpz_vec_sum_max_bits(&sum_bits, &max_bits, a, len);

        fmpz_zero(max);
        fmpz_zero(sum);
        for (j = 0; j < len; j++)
        {
            fmpz_abs(t, a + j);
            fmpz_add(sum, sum, t);
            if (fmpz_cmp(max, t) < 0)
                fmpz_set(max, t);
        }

        if (sum_bits != fmpz_bits(sum))
        {
            flint_printf("FAIL: sum bits is wrong\n");
            fflush(stdout);
            flint_abort();
        }

        if (max_bits != fmpz_bits(max))
        {
            flint_printf("FAIL: max bits is wrong\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        fmpz_clear(max);
        fmpz_clear(sum);
        fmpz_clear(t);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
