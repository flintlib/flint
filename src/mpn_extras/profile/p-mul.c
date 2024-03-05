/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "mpn_extras.h"
#include "profiler.h"

#define N_MAX 16
#define N_MAX2 256

static void measure(mp_ptr rpg, mp_ptr rpf, mp_ptr xp, mp_ptr yp, slong m)
{
    double t1, t2, FLINT_SET_BUT_UNUSED(__);
    slong n;

    flint_printf("m = %3wd:", m);

    for (n = 1; n <= m; n = (n <= N_MAX ? n + 1 : n * 1.2))
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

int main(int argc, char ** argv)
{
    mp_limb_t xp[N_MAX2], yp[N_MAX2], rpf[2 * N_MAX2], rpg[2 * N_MAX2];
    slong m;

    if (argc == 2)
    {
        measure(rpg, rpf, xp, yp, strtol(argv[1], NULL, 10));
        goto end;
    }

    flint_printf("mpn_mul vs flint_mpn_mul\n\n");

    for (m = 1; m <= N_MAX2; m = (m <= N_MAX ? m + 1 : m * 1.2))
        measure(rpg, rpf, xp, yp, m);

end:
    flint_cleanup_master();
    return 0;
}
