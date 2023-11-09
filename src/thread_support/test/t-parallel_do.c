/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "fmpz.h"

typedef struct
{
    int * res;
}
f_param_t;

void
f(slong i, void * param)
{
    f_param_t * p = (f_param_t *) param;

    p->res[i] = i * i;
}

TEST_FUNCTION_START(thread_support_parallel_do, state)
{
    slong iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        int * resx;
        int * resy;
        slong i, n;
        f_param_t workx, worky;

        n = n_randint(state, 1000);

        flint_set_num_threads(n_randint(state, 10) + 1);

        resx = flint_malloc(n * sizeof(int));
        resy = flint_malloc(n * sizeof(int));

        workx.res = resx;
        worky.res = resy;

        flint_parallel_do(f, &workx, n, n_randint(state, 5), FLINT_PARALLEL_UNIFORM);
        flint_parallel_do(f, &worky, n, n_randint(state, 5), FLINT_PARALLEL_STRIDED);

        for (i = 0; i < n; i++)
        {
            if (resx[i] != resy[i] || resx[i] != i * i)
            {
                flint_printf("FAIL\n");
                flint_printf("num_threads = %wd, i = %wd/%wd\n", flint_get_num_threads(), i, n);
                flint_abort();
            }
        }

        flint_free(resx);
        flint_free(resy);
    }

    TEST_FUNCTION_END(state);
}
