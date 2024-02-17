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

#if !FLINT_HAVE_ADX
# error
#endif

#define N_MIN 1
#define N_MAX 128

#define mpn_sqr_basecase __gmpn_sqr_basecase
void mpn_sqr_basecase(mp_ptr, mp_srcptr, mp_size_t);

int main(void)
{
    mp_limb_t rf[N_MAX];
    mp_limb_t rg[2 * N_MAX];
    mp_limb_t xp[N_MAX];
    mp_size_t n;

    flint_printf("mpn_sqr_basecase / flint_mpn_sqrhigh_basecase:\n");
    for (n = N_MIN; n <= N_MAX; n++)
    {
        double t1, t2, __attribute__((unused)) __;
        flint_printf("n = %2wd:", n);

        mpn_random2(xp, n);

        TIMEIT_START
        flint_mpn_sqrhigh_basecase(rf, xp, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        flint_mpn_sqr(rg, xp, n);
        TIMEIT_STOP_VALUES(__, t2)

        flint_printf("%7.2fx\n", t2 / t1);
    }

    flint_cleanup_master();

    return 0;
}

#undef N_MAX
#undef N_MIN
