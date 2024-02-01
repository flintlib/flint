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

#define mpn_mul_basecase __gmpn_mul_basecase
void mpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
void mpfr_mulhigh_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

#define N_MAX FLINT_MPN_MULHIGH_N_FUNC_TAB_WIDTH

int main(void)
{
    mp_limb_t rf[N_MAX];
    mp_limb_t rg[2 * N_MAX];
    mp_limb_t rm[2 * N_MAX];
    mp_limb_t xp[N_MAX];
    mp_limb_t yp[N_MAX];
    mp_size_t n;

    flint_printf("           GMP       MPFR\n");
    for (n = 1; n <= N_MAX; n++)
    {
        double t1, t2, t3, __attribute__((unused)) __;
        flint_printf("n = %2wd:", n);

        mpn_random2(xp, n);
        mpn_random2(yp, n);

        TIMEIT_START
        flint_mpn_mulhigh_n(rf, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        mpn_mul_basecase(rg, xp, n, yp, n);
        TIMEIT_STOP_VALUES(__, t2)

        TIMEIT_START
        mpfr_mulhigh_n(rm, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t3)

        flint_printf("%7.2fx  %7.2fx\n", t2 / t1, t3 / t1);
    }

    flint_cleanup_master();

    return 0;
}
