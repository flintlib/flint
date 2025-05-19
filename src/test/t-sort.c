/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "flint.h"
#include "perm.h"

static int
cmp_first(const void * _i, const void * _j, void * _imax)
{
    const slong * i = _i, * j = _j;
    slong * imax = _imax;
    *imax = FLINT_MAX(*imax, FLINT_MAX(*i, *j));
    return *i - *j;
}

TEST_FUNCTION_START(flint_sort, state)
{
    for (slong i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong imax = -1, imax1 = -1;

        slong n = n_randint(state, 10);

        slong * vec = _perm_init(2*n);
        _perm_randtest(vec, 2*n, state);

        int dup = n_randint(state, 2);
        if (dup && n > 0)
        {
            slong m = n_randint(state, n);
            _perm_randtest(vec + 2*m, 2*n - 2*m, state);

            for (slong k = 0; k < n; k++)
                vec[2*k + 1] = k;
        }

        slong * vec1 = flint_malloc(2*n*sizeof(slong));
        memcpy(vec1, vec, 2*n*sizeof(slong));

        flint_sort(vec, n, 2*sizeof(slong), cmp_first, &imax);
        flint_merge_sort(vec1, n, 2*sizeof(slong), cmp_first, &imax1);

        if (n > 0 && memcmp(vec, vec1, 2*sizeof(slong)))
            TEST_FUNCTION_FAIL("inconsistent results\n");

        if (imax != imax1)
            TEST_FUNCTION_FAIL("inconsistent data\n");

        for (slong k = 0; k < n - 1; k++)
        {
            if (vec[2 * k] > vec[2 * k + 2])
                TEST_FUNCTION_FAIL("wrong order\n");
            if (dup && vec[2 * k] == vec[2 * k + 2]
                    && vec[2 * k + 1] >= vec[2 * k + 3])
                TEST_FUNCTION_FAIL("stability\n");
        }

        if (!dup)
        {
            flint_merge_sort(vec, 2 * n, sizeof(slong), cmp_first, &imax);

            for (slong k = 0; k < 2 * n; k++)
            {
                if (!dup && vec[k] != k)
                    TEST_FUNCTION_FAIL("missing entry\n");
            }

            if (imax != 2 * n - 1)
                TEST_FUNCTION_FAIL("data\n");
        }

        flint_free(vec1);
        _perm_clear(vec);
    }

    TEST_FUNCTION_END(state);
}

