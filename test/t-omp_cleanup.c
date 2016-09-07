/*
    Authored 2016 by Daniel S. Roche

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

int thread_active = 0;
#pragma omp threadprivate(thread_active)

int thread_count = 0;
int free_count = 0;

void cleanup()
{
    if (thread_active)
    {
#pragma omp atomic
        --thread_count;
#pragma omp atomic
        ++free_count;
        thread_active = 0;
    }
}

void foo()
{
    if (!thread_active)
    {
#pragma omp atomic
        ++thread_count;
        flint_register_cleanup_function(cleanup);
        thread_active = 1;
    }
}

int main(void)
{
#ifdef _OPENMP
    int i, result;
#endif
    FLINT_TEST_INIT(state);

    flint_printf("omp_cleanup....");
    fflush(stdout);

#ifdef _OPENMP

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        int threads, tmin;
        slong iters;

        threads = (int)n_randint(state, 8) + 1;
        flint_set_num_threads(threads);

        iters = (slong)n_randint(state, 1000) + threads;

        free_count = 0;

#pragma omp parallel
        {
            slong j;

#pragma omp for
            for (j = 0; j < iters; ++j)
                foo();

            flint_parallel_cleanup();
        }

#ifdef _OPENMP
        tmin = FLINT_MIN(1, threads - 1);
#else
        tmin = 0;
#endif

        result = (thread_count == 1) && (free_count >= tmin);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("tcount = %d, fcount = %d\n",
                thread_count, free_count);
            abort();
        }
    }

    FLINT_TEST_CLEANUP(state);

    if (thread_count != 0)
    {
        flint_printf("FAIL: final cleanup\n");
        abort();
    }
    
    flint_printf("PASS\n");

#else
    FLINT_TEST_CLEANUP(state);
    flint_printf("SKIPPED\n");
#endif

    return 0;
}
