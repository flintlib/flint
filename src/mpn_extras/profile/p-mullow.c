/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "profiler.h"

void __gmpn_mullo_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

#define N_MAX 2100

int
main(void)
{
    mp_limb_t rp[2 * N_MAX];
    mp_limb_t xp[N_MAX];
    mp_limb_t yp[N_MAX];
    mp_size_t n;

    flint_printf("        speedup\n");
    flint_printf("   n    vs mpn_mullo_n     vs mul    vs basecase\n");
    for (n = 1; n < N_MAX; n = FLINT_MAX(n + 1, n * 1.05))
    {
        double t1, t2, t3, t4, FLINT_SET_BUT_UNUSED(__);
        flint_printf("%4wd", n);

        mpn_random2(xp, n);
        mpn_random2(yp, n);

        TIMEIT_START
        __gmpn_mullo_n(rp, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        flint_mpn_mul_n(rp, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t2)

        if (n >= 9)
        {
            TIMEIT_START
            flint_mpn_mullow_basecase(rp, xp, yp, n);
            TIMEIT_STOP_VALUES(__, t3)
        }
        else
            t3 = 0.0;

        TIMEIT_START
        flint_mpn_mullow_n(rp, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t4)

        flint_printf("      %7.2fx       %7.2f       %7.2f\n", t1 / t4, t2 / t4, t3 / t4);
    }

    flint_cleanup_master();

    return 0;
}

#undef N_MAX