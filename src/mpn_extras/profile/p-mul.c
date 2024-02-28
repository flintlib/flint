/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "profiler.h"

#define MAXN 16
#define MAXN2 256

int main(void)
{
    slong m, n;

    flint_printf("mpn_mul vs flint_mpn_mul\n\n");

    for (m = 1; m <= MAXN2; m = (m <= MAXN ? m + 1 : m * 1.2))
    {
        mp_limb_t xp[MAXN2], yp[MAXN2], rpf[2 * MAXN2], rpg[2 * MAXN2];
        double t1, t2, FLINT_SET_BUT_UNUSED(__);

        flint_printf("m = %3wd:", n);

        for (n = 1; n <= m; n = (n <= MAXN ? n + 1 : n * 1.2))
        {
            mpn_random2(xp, m);
            mpn_random2(yp, n);

            TIMEIT_START
            mpn_mul(rpg, xp, m, yp, n);
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            flint_mpn_mul(rpf, xp, m, yp, n);
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf(" %.2f", t1 / t2);
        }

        flint_printf("\n");
    }

    flint_cleanup_master();
    return 0;
}
