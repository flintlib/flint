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

#define N_MAX FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH

int main(void)
{
    mp_limb_t rn[N_MAX];
    mp_limb_t ru[N_MAX];
    mp_limb_t xp[N_MAX];
    mp_limb_t yp[N_MAX];
    mp_size_t n;

    flint_printf("Times: Normalised / Unnormalised\n");
    for (n = 1; n <= N_MAX; n++)
    {
        flint_printf("n = %2wd:", n);

        for (slong ix = 0; ix < 10; ix++)
        {
            double t1, t2, __attribute__((unused)) __;

            mpn_random2(xp, n);
            mpn_random2(yp, n);
            xp[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));
            yp[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));

            TIMEIT_START
            flint_mpn_mulhigh(ru, xp, yp, n);
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            flint_mpn_mulhigh_normalised(rn, xp, yp, n);
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%7.2fx", t2 / t1);
        }
        flint_printf("\n");
    }

    flint_cleanup_master();

    return 0;
}
