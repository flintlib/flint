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

static void measure(mp_ptr rpg, mp_ptr rpf, mp_ptr xp, mp_ptr yp, slong m, slong n_max)
{
    double t1, t2, FLINT_SET_BUT_UNUSED(__);
    slong n;

    flint_printf("m = %3wd:", m);

    for (n = 1; n <= n_max; n = (n <= N_MAX ? n + 1 : n * 1.2))
    {
        mp_ptr ap, bp;
        slong mt, nt;

        if (m >= n)
        {
            mt = m;
            nt = n;
            ap = xp;
            bp = yp;
        }
        else
        {
            mt = n;
            nt = m;
            ap = yp;
            bp = xp;
        }

        mpn_random2(ap, mt);
        mpn_random2(bp, nt);

        TIMEIT_START
        mpn_mul(rpg, ap, mt, bp, nt);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        flint_mpn_mul(rpf, ap, mt, bp, nt);
        TIMEIT_STOP_VALUES(__, t2)

        flint_printf(" %.2f", t1 / t2);
    }

    flint_printf("\n");
}

int main(int argc, char ** argv)
{
    mp_ptr xp, yp, rpf, rpg;
    slong m;

    xp = flint_malloc(sizeof(mp_limb_t) * N_MAX2);
    yp = flint_malloc(sizeof(mp_limb_t) * N_MAX2);
    rpf = flint_malloc(2 * sizeof(mp_limb_t) * N_MAX2);
    rpg = flint_malloc(2 * sizeof(mp_limb_t) * N_MAX2);

    if (argc == 2)
    {
        measure(rpg, rpf, xp, yp, strtol(argv[1], NULL, 10), N_MAX2);
        goto end;
    }
    else if (argc == 3)
    {
        slong m = strtol(argv[1], NULL, 10);
        slong n_max = strtol(argv[2], NULL, 10);
        xp = flint_realloc(xp, sizeof(mp_limb_t) * m);
        yp = flint_realloc(yp, sizeof(mp_limb_t) * n_max);
        rpf = flint_realloc(rpf, sizeof(mp_limb_t) * (m + n_max));
        rpg = flint_realloc(rpg, sizeof(mp_limb_t) * (m + n_max));
        measure(rpg, rpf, xp, yp, m, n_max);
        goto end;
    }

    flint_printf("mpn_mul vs flint_mpn_mul\n\n");

    for (m = 1; m <= N_MAX2; m = (m <= N_MAX ? m + 1 : m * 1.2))
        measure(rpg, rpf, xp, yp, m, m);

end:
    flint_free(xp);
    flint_free(yp);
    flint_free(rpf);
    flint_free(rpg);

    flint_cleanup_master();
    return 0;
}
